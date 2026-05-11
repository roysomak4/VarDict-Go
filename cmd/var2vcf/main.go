package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

const Version = "1.8.2"

type Config struct {
	SampleName    string
	RefPath       string
	BedPath       string
	MinTotalDepth int
	MinVarDepth   int
	MinAF         float64
	MinPMean      float64
	MinQual       float64
	MinMapQ       float64
	MinSN         float64
	GTFreq        float64
	MaxMSI        int
	MaxNM         float64
	ClusterBP     int
	FilterPStd    int  // 0=off, 1=on; default 1
	MinSplitReads int
	NoEnd         bool // -E flag: if set, omit END tag from INFO
	AllVariants   bool
	Strict        bool
}

// VariantRecord maps the 40-column post-teststrandbias VarDict output.
// Columns 20-21 are the Fisher p-value and odds ratio inserted by teststrandbias;
// everything from col 22 onward is the original VarDict col 20+ shifted right by 2.
type VariantRecord struct {
	Sample   string
	Gene     string
	Chr      string
	Start    int
	End      int
	Ref      string
	Alt      string
	DP       int     // col 7: total depth
	VD       int     // col 8: variant depth (raw read count)
	RFwd     int     // col 9
	RRev     int     // col 10
	VFwd     int     // col 11
	VRev     int     // col 12
	Genotype string  // col 13
	AF       float64 // col 14
	Bias     string  // col 15
	PMean    float64 // col 16
	PStd     float64 // col 17
	Qual     float64 // col 18
	QStd     float64 // col 19
	SBF      float64 // col 20: Fisher p-value (inserted by teststrandbias)
	OddRatio float64 // col 21: odds ratio (inserted by teststrandbias)
	MapQ     float64 // col 22
	SN       float64 // col 23
	HiAF     float64 // col 24
	AdjAF    float64 // col 25
	Shift3   int     // col 26
	MSI      float64 // col 27
	MSILen   float64 // col 28
	NM       float64 // col 29
	HiCnt    int     // col 30: high-quality variant reads (used for depth filter)
	HiCov    int     // col 31: high-quality total reads
	LSeq     string  // col 32: 5' flanking sequence
	RSeq     string  // col 33: 3' flanking sequence
	// col 34 (seg) is not used in output
	Type    string // col 35
	GAmp    string // col 36: DUPRATE for non-amplicon; GDAMP for amplicon
	TAmp    string // col 37: SV split/span info or total amplicons
	NCAmp   string // col 38: CRISPR for non-amplicon; NCAMP for amplicon
	AmpFlag string // col 39: amplicon mismatch flag; presence signals amplicon mode
	IsAmp   bool   // true when col 39 is present (len(fields) > 39)
}

func main() {
	config := parseFlags()

	variants, sampleName, err := readAndGroupVariants(os.Stdin)
	if err != nil {
		log.Fatalf("Error reading variants: %v", err)
	}

	if config.SampleName != "" {
		sampleName = config.SampleName
	}
	if sampleName == "" {
		sampleName = "SAMPLE"
	}

	printVCFHeader(config, sampleName, os.Stdout)
	processVariants(variants, config, os.Stdout)
}

func parseFlags() Config {
	c := Config{}
	flag.StringVar(&c.SampleName, "N", "", "Sample name override")
	flag.StringVar(&c.RefPath, "G", "", "Path to reference FASTA (for ##reference header)")
	flag.StringVar(&c.BedPath, "b", "", "Path to BED file (for ##contig headers)")
	flag.IntVar(&c.MinTotalDepth, "d", 3, "Minimum total depth")
	flag.IntVar(&c.MinVarDepth, "v", 2, "Minimum high-quality variant depth")
	flag.Float64Var(&c.MinAF, "f", 0.02, "Minimum allele frequency")
	flag.Float64Var(&c.MinPMean, "p", 8.0, "Minimum mean position in read")
	flag.Float64Var(&c.MinQual, "q", 22.5, "Minimum mean base quality")
	flag.Float64Var(&c.MinMapQ, "Q", 10.0, "Minimum mean mapping quality")
	flag.Float64Var(&c.MinSN, "o", 1.5, "Minimum signal to noise ratio")
	flag.Float64Var(&c.GTFreq, "F", 0.2, "Genotype frequency threshold for homozygous call")
	flag.IntVar(&c.MaxMSI, "I", 12, "Maximum non-monomer MSI for AF<0.5 variants")
	flag.Float64Var(&c.MaxNM, "m", 5.25, "Maximum mean mismatches in reads")
	flag.IntVar(&c.ClusterBP, "c", 0, "Filter SNVs within this many bp of each other (0=disabled)")
	flag.IntVar(&c.FilterPStd, "P", 1, "Filter variants with pstd=0: 1=yes (default), 0=no")
	flag.IntVar(&c.MinSplitReads, "T", 1, "Minimum split reads for structural variants")
	flag.BoolVar(&c.NoEnd, "E", false, "If set, do not print END tag in INFO field")
	flag.BoolVar(&c.AllVariants, "A", false, "Output all variants at the same position")
	flag.BoolVar(&c.Strict, "S", false, "Only output variants that pass all filters")
	flag.Parse()
	return c
}

func readAndGroupVariants(reader io.Reader) (map[string]map[int][]*VariantRecord, string, error) {
	variants := make(map[string]map[int][]*VariantRecord)
	var sampleName string

	scanner := bufio.NewScanner(reader)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.Contains(line, "R_HOME") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 36 {
			continue
		}
		record, err := parseVariantRecord(fields)
		if err != nil {
			continue
		}
		if sampleName == "" {
			sampleName = record.Sample
		}
		if variants[record.Chr] == nil {
			variants[record.Chr] = make(map[int][]*VariantRecord)
		}
		variants[record.Chr][record.Start] = append(variants[record.Chr][record.Start], record)
	}
	return variants, sampleName, scanner.Err()
}

func parseVariantRecord(fields []string) (*VariantRecord, error) {
	atoi := func(s string) int { v, _ := strconv.Atoi(s); return v }
	atof := func(s string) float64 { v, _ := strconv.ParseFloat(s, 64); return v }
	get := func(i int) string {
		if i < len(fields) {
			return fields[i]
		}
		return ""
	}

	r := &VariantRecord{
		Sample:   fields[0],
		Gene:     fields[1],
		Chr:      fields[2],
		Start:    atoi(fields[3]),
		End:      atoi(fields[4]),
		Ref:      fields[5],
		Alt:      fields[6],
		DP:       atoi(fields[7]),
		VD:       atoi(fields[8]),
		RFwd:     atoi(fields[9]),
		RRev:     atoi(fields[10]),
		VFwd:     atoi(fields[11]),
		VRev:     atoi(fields[12]),
		Genotype: fields[13],
		AF:       atof(fields[14]),
		Bias:     fields[15],
		PMean:    atof(fields[16]),
		PStd:     atof(fields[17]),
		Qual:     atof(fields[18]),
		QStd:     atof(fields[19]),
		SBF:      atof(fields[20]),
		OddRatio: atof(fields[21]),
		MapQ:     atof(fields[22]),
		SN:       atof(fields[23]),
		HiAF:     atof(fields[24]),
		AdjAF:    atof(fields[25]),
		Shift3:   atoi(get(26)),
		MSI:      atof(get(27)),
		MSILen:   atof(get(28)),
		NM:       atof(get(29)),
		HiCnt:    atoi(get(30)),
		HiCov:    atoi(get(31)),
		LSeq:     get(32),
		RSeq:     get(33),
		// col 34 (seg) skipped
		Type:    get(35),
		GAmp:    get(36),
		TAmp:    get(37),
		NCAmp:   get(38),
		AmpFlag: get(39),
		IsAmp:   len(fields) > 39,
	}

	// Normalize odds ratio once at parse time: Inf→0, (0,1)→invert
	if math.IsInf(r.OddRatio, 1) || math.IsInf(r.OddRatio, -1) {
		r.OddRatio = 0
	} else if r.OddRatio > 0 && r.OddRatio < 1 {
		r.OddRatio = 1 / r.OddRatio
	}

	return r, nil
}

func processVariants(variants map[string]map[int][]*VariantRecord, config Config, writer io.Writer) {
	for _, chr := range sortChromosomes(variants) {
		processChromosome(chr, variants[chr], config, writer)
	}
}

func processChromosome(chr string, positions map[int][]*VariantRecord, config Config, writer io.Writer) {
	posList := make([]int, 0, len(positions))
	for pos := range positions {
		posList = append(posList, pos)
	}
	sort.Ints(posList)

	// Deferred printing: buffer the current record so the Cluster filter can
	// retroactively update the previous record's filter before it is printed.
	var pendingCols1, pendingFilter, pendingCols2 string

	flush := func() {
		if pendingCols1 == "" {
			return
		}
		if !config.Strict || pendingFilter == "PASS" {
			fmt.Fprintf(writer, "%s\t%s\t%s\n", pendingCols1, pendingFilter, pendingCols2)
		}
		pendingCols1 = ""
	}

	var prevSNVStart int

	for _, pos := range posList {
		records := positions[pos]

		// Sort by AF descending (highest AF first, matching Perl behaviour)
		sort.Slice(records, func(i, j int) bool {
			return records[i].AF > records[j].AF
		})

		limit := 1
		if config.AllVariants {
			limit = len(records)
		}

		seen := make(map[string]bool)

		for i := 0; i < limit; i++ {
			r := records[i]
			if r.Ref == "" {
				continue
			}
			key := fmt.Sprintf("%s-%d-%d-%s-%s", r.Chr, r.Start, r.End, r.Ref, r.Alt)
			if seen[key] {
				continue
			}
			seen[key] = true

			if r.Type == "" {
				r.Type = "REF"
			}

			filters, clusterApplied := applyFilters(r, config, prevSNVStart)

			// SV without enough split-read support: silently skip and clear buffer
			if len(filters) == 1 && filters[0] == "SKIP" {
				pendingCols1 = ""
				continue
			}

			filterStr := "PASS"
			if len(filters) > 0 {
				filterStr = strings.Join(filters, ";")
			}

			// Strict mode: don't process non-PASS records at all
			if config.Strict && filterStr != "PASS" {
				continue
			}

			// Retroactively mark the buffered previous SNV with the Cluster filter
			if clusterApplied && pendingCols1 != "" {
				pendingFilter = fmt.Sprintf("Cluster%dbp", config.ClusterBP)
			}

			flush()

			cols1, cols2 := buildVCFRecordParts(r, config)
			pendingCols1 = cols1
			pendingFilter = filterStr
			pendingCols2 = cols2

			if r.Type == "SNV" && filterStr == "PASS" {
				prevSNVStart = r.Start
			}
		}
	}

	flush()
}

// applyFilters returns the list of filter tags for the record.
// Returns ["SKIP"] when the record should be silently dropped (SV without support).
// clusterApplied is true when the Cluster filter was added, so the caller can
// retroactively apply it to the previously buffered record as well.
func applyFilters(r *VariantRecord, config Config, prevSNVStart int) ([]string, bool) {
	var filters []string

	// Depth filters use HiCnt (high-quality variant reads, col 30), not raw VD
	if r.DP < config.MinTotalDepth {
		if !(float64(r.HiCnt)*r.HiAF >= 0.5) {
			filters = append(filters, fmt.Sprintf("d%d", config.MinTotalDepth))
		}
	}
	if r.HiCnt < config.MinVarDepth {
		if !(float64(r.HiCnt)*r.HiAF >= 0.5) {
			filters = append(filters, fmt.Sprintf("v%d", config.MinVarDepth))
		}
	}

	if r.AF < config.MinAF {
		filters = append(filters, fmt.Sprintf("f%g", config.MinAF))
	}
	if r.PMean < config.MinPMean {
		filters = append(filters, fmt.Sprintf("p%g", config.MinPMean))
	}
	if config.FilterPStd == 1 && r.PStd == 0 && !r.IsAmp && r.AF < 0.35 {
		filters = append(filters, "pSTD")
	}
	if r.Qual < config.MinQual {
		filters = append(filters, fmt.Sprintf("q%g", config.MinQual))
	}
	if r.MapQ < config.MinMapQ && r.AF < 0.8 {
		filters = append(filters, fmt.Sprintf("Q%g", config.MinMapQ))
	}
	if r.SN < config.MinSN {
		filters = append(filters, fmt.Sprintf("SN%g", config.MinSN))
	}
	if r.NM > config.MaxNM {
		filters = append(filters, fmt.Sprintf("NM%g", config.MaxNM))
	}

	// MSI filter
	lenDiff := absInt(len(r.Ref) - len(r.Alt))
	if (r.MSI > float64(config.MaxMSI) && r.MSILen > 1 && r.AF < 0.2 && float64(lenDiff) == r.MSILen) ||
		(r.MSI >= 13 && r.MSILen == 1 && r.AF <= 0.275 && float64(lenDiff) == r.MSILen) {
		filters = append(filters, fmt.Sprintf("MSI%d", config.MaxMSI))
	}

	// Strand bias filter — matches Perl's exact 5-condition check
	if r.HiAF < 0.25 && r.Bias == "2;1" && r.SBF < 0.01 &&
		(r.OddRatio > 5 || r.OddRatio == 0) && (r.End-r.Start) < 100 {
		filters = append(filters, "Bias")
	}

	// LongMSI filter
	if float64(lenDiff) == r.MSILen {
		if r.HiAF <= 0.275 && r.MSI >= 13 {
			filters = append(filters, "LongMSI")
		} else if r.HiAF <= 0.2 && r.MSI >= 8 && r.MSILen > 1 {
			filters = append(filters, "LongMSI")
		}
	}

	// Cluster filter: only when -c > 0 and no other filters triggered
	clusterApplied := false
	if config.ClusterBP > 0 && r.Type == "SNV" && len(filters) == 0 &&
		prevSNVStart > 0 && r.Start-prevSNVStart < config.ClusterBP {
		filters = append(filters, fmt.Sprintf("Cluster%dbp", config.ClusterBP))
		clusterApplied = true
	}

	// SV without sufficient split-read support: signal caller to skip silently
	if strings.Contains(r.Alt, "<") && !r.IsAmp {
		sr, _, _ := parseTAmp(r.TAmp)
		if sr < config.MinSplitReads {
			return []string{"SKIP"}, false
		}
	}

	return filters, clusterApplied
}

// buildVCFRecordParts returns the two tab-separated column groups that bracket
// the FILTER field: (CHROM POS ID REF ALT QUAL) and (INFO FORMAT SAMPLE).
// The caller assembles: cols1 + "\t" + filter + "\t" + cols2 + "\n".
func buildVCFRecordParts(r *VariantRecord, config Config) (cols1, cols2 string) {
	ref := r.Ref
	alt := r.Alt
	gt := calculateGenotype(r, config)

	if ref == alt {
		alt = "."
		gt = "0/0"
	}

	// QUAL: 0 when vd<=1, otherwise floor(log2(vd) * qual)
	var vcfQual int
	if r.VD > 1 {
		vcfQual = int(math.Log2(float64(r.VD)) * r.Qual)
	}

	cols1 = fmt.Sprintf("%s\t%d\t.\t%s\t%s\t%d", r.Chr, r.Start, ref, alt, vcfQual)
	cols2 = fmt.Sprintf("%s\t%s\t%s", buildINFO(r, config), "GT:DP:VD:AD:AF:RD:ALD", buildSampleData(r, gt))
	return
}

func calculateGenotype(r *VariantRecord, config Config) string {
	if r.Ref == r.Alt {
		return "0/0"
	}
	if 1-r.AF < config.GTFreq {
		return "1/1"
	}
	if r.AF >= 0.5 {
		return "1/0"
	}
	if r.AF >= config.MinAF {
		return "0/1"
	}
	return "0/0"
}

func buildINFO(r *VariantRecord, config Config) string {
	// Bias string: replace ';' with ':' to avoid breaking semicolon-delimited INFO
	bias := strings.ReplaceAll(r.Bias, ";", ":")
	sampleName := strings.ReplaceAll(r.Sample, " ", "_")

	var b strings.Builder
	b.WriteString(fmt.Sprintf("SAMPLE=%s", sampleName))
	b.WriteString(fmt.Sprintf(";TYPE=%s", r.Type))
	b.WriteString(fmt.Sprintf(";DP=%d", r.DP))
	if !config.NoEnd {
		b.WriteString(fmt.Sprintf(";END=%d", r.End))
	}
	b.WriteString(fmt.Sprintf(";VD=%d", r.VD))
	b.WriteString(fmt.Sprintf(";AF=%g", r.AF))
	b.WriteString(fmt.Sprintf(";BIAS=%s", bias))
	b.WriteString(fmt.Sprintf(";REFBIAS=%d:%d", r.RFwd, r.RRev))
	b.WriteString(fmt.Sprintf(";VARBIAS=%d:%d", r.VFwd, r.VRev))
	b.WriteString(fmt.Sprintf(";PMEAN=%g", r.PMean))
	b.WriteString(fmt.Sprintf(";PSTD=%g", r.PStd))
	b.WriteString(fmt.Sprintf(";QUAL=%g", r.Qual))
	b.WriteString(fmt.Sprintf(";QSTD=%g", r.QStd))
	b.WriteString(fmt.Sprintf(";SBF=%g", r.SBF))
	b.WriteString(fmt.Sprintf(";ODDRATIO=%g", r.OddRatio))
	b.WriteString(fmt.Sprintf(";MQ=%g", r.MapQ))
	b.WriteString(fmt.Sprintf(";SN=%g", r.SN))
	b.WriteString(fmt.Sprintf(";HIAF=%g", r.HiAF))
	b.WriteString(fmt.Sprintf(";ADJAF=%g", r.AdjAF))
	b.WriteString(fmt.Sprintf(";SHIFT3=%d", r.Shift3))
	b.WriteString(fmt.Sprintf(";MSI=%g", r.MSI))
	b.WriteString(fmt.Sprintf(";MSILEN=%g", r.MSILen))
	b.WriteString(fmt.Sprintf(";NM=%g", r.NM))
	b.WriteString(fmt.Sprintf(";HICNT=%d", r.HiCnt))
	b.WriteString(fmt.Sprintf(";HICOV=%d", r.HiCov))
	b.WriteString(fmt.Sprintf(";LSEQ=%s", r.LSeq))
	b.WriteString(fmt.Sprintf(";RSEQ=%s", r.RSeq))

	// Amplicon mode fields (TODO: full amplicon mode support)
	if r.IsAmp {
		b.WriteString(fmt.Sprintf(";GDAMP=%s;TLAMP=%s;NCAMP=%s;AMPFLAG=%s",
			r.GAmp, r.TAmp, r.NCAmp, r.AmpFlag))
	} else {
		if r.GAmp != "" {
			b.WriteString(fmt.Sprintf(";DUPRATE=%s", r.GAmp))
		}
		if r.NCAmp != "" {
			b.WriteString(fmt.Sprintf(";CRISPR=%s", r.NCAmp))
		}
	}

	// SV fields (non-amplicon only)
	if !r.IsAmp {
		sr, sp, _ := parseTAmp(r.TAmp)
		if strings.Contains(r.Alt, "<") {
			svType := strings.Trim(r.Alt, "<>")
			svLen := r.End - r.Start
			if strings.Contains(svType, "INV") {
				svLen++
			}
			b.WriteString(fmt.Sprintf(";SVTYPE=%s;SVLEN=%d", svType, svLen))
		}
		if r.TAmp != "" {
			b.WriteString(fmt.Sprintf(";SPLITREAD=%d;SPANPAIR=%d", sr, sp))
		}
	}

	return b.String()
}

func buildSampleData(r *VariantRecord, gt string) string {
	rd := r.RFwd + r.RRev

	// AD: reference depth only when VD=0, otherwise ref,var
	var ad string
	if r.VD == 0 {
		ad = strconv.Itoa(rd)
	} else {
		ad = fmt.Sprintf("%d,%d", rd, r.VD)
	}

	return fmt.Sprintf("%s:%d:%d:%s:%g:%d,%d:%d,%d",
		gt,
		r.DP,
		r.VD,
		ad,
		r.AF,
		r.RFwd, r.RRev,
		r.VFwd, r.VRev,
	)
}

// parseTAmp extracts splitReads, spanPairs, cluster from a "-"-delimited TAmp field.
func parseTAmp(tamp string) (splitReads, spanPairs, cluster int) {
	if !strings.Contains(tamp, "-") {
		return
	}
	parts := strings.SplitN(tamp, "-", 3)
	splitReads, _ = strconv.Atoi(parts[0])
	if len(parts) > 1 {
		spanPairs, _ = strconv.Atoi(parts[1])
	}
	if len(parts) > 2 {
		cluster, _ = strconv.Atoi(parts[2])
	}
	return
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

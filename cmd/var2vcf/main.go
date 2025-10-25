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

const Version = "1.0.0"

// Config holds command-line configuration
type Config struct {
	SampleName      string
	RefPath         string
	BedPath         string
	MinTotalDepth   int
	MinVariantDepth int
	MinAF           float64
	MinPMean        float64
	MinQual         float64
	MinMapQ         float64
	MinSN           float64
	GTFreq          float64
	PrintEndTag     bool
	AllVariants     bool
	FilterPStd      bool
	MinSplitReads   int
}

// VariantRecord represents a single variant from the input
type VariantRecord struct {
	Sample    string
	Gene      string
	Chr       string
	Start     int
	End       int
	Ref       string
	Alt       string
	Depth     int
	AltDepth  int
	RefFwd    int
	RefRev    int
	AltFwd    int
	AltRev    int
	Genotype  string
	AF        float64
	Bias      string
	PMean     float64
	PStd      float64
	Qual      float64
	QStd      float64
	MapQ      float64
	QRatio    float64
	HiAF      float64
	AdjAF     float64
	NM        float64
	MQ        float64
	DupRate   float64
	SV        string
	Type      string
	Status    string
	PValue    float64
	OddsRatio float64
}

func main() {
	config := parseFlags()

	// Read and group variants by chromosome
	variants, sampleName, err := readAndGroupVariants(os.Stdin)
	if err != nil {
		log.Fatalf("Error reading variants: %v", err)
	}

	// Override sample name if provided via -N flag
	if config.SampleName != "" {
		sampleName = config.SampleName
	}

	if sampleName == "" {
		sampleName = "SAMPLE"
	}

	// Print VCF header
	printVCFHeader(config, sampleName, os.Stdout)

	// Process variants chromosome by chromosome
	processVariants(variants, config, os.Stdout)
}

func parseFlags() Config {
	config := Config{}

	flag.StringVar(&config.SampleName, "N", "", "Sample name")
	flag.StringVar(&config.RefPath, "G", "", "Reference genome path")
	flag.StringVar(&config.BedPath, "b", "", "BED file path for contig information")
	flag.IntVar(&config.MinTotalDepth, "d", 8, "Minimum total depth")
	flag.IntVar(&config.MinVariantDepth, "v", 4, "Minimum variant depth")
	flag.Float64Var(&config.MinAF, "f", 0.01, "Minimum allele frequency")
	flag.Float64Var(&config.MinPMean, "p", 5.0, "Minimum mean position in read")
	flag.Float64Var(&config.MinQual, "q", 22.5, "Minimum base quality")
	flag.Float64Var(&config.MinMapQ, "Q", 20.0, "Minimum mapping quality")
	flag.Float64Var(&config.MinSN, "S", 1.5, "Minimum signal-to-noise ratio")
	flag.Float64Var(&config.GTFreq, "F", 0.2, "Genotype frequency threshold")
	flag.BoolVar(&config.PrintEndTag, "E", false, "Print END tag in INFO field")
	flag.BoolVar(&config.AllVariants, "A", false, "Print all variants at the same position")
	flag.BoolVar(&config.FilterPStd, "P", false, "Filter variants with position std dev = 0")
	flag.IntVar(&config.MinSplitReads, "T", 1, "Minimum split reads for structural variants")

	flag.Parse()

	return config
}

func readAndGroupVariants(reader io.Reader) (map[string]map[int][]*VariantRecord, string, error) {
	variants := make(map[string]map[int][]*VariantRecord)
	var sampleName string
	scanner := bufio.NewScanner(reader)

	// Increase buffer size for long lines
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)

	for scanner.Scan() {
		line := scanner.Text()

		// Skip R environment lines
		if strings.Contains(line, "R_HOME") {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 32 {
			continue // Skip malformed lines
		}

		record, err := parseVariantRecord(fields)
		if err != nil {
			continue // Skip unparseable lines
		}

		// Store first sample name
		if sampleName == "" {
			sampleName = record.Sample
		}

		// Group by chr -> pos
		if variants[record.Chr] == nil {
			variants[record.Chr] = make(map[int][]*VariantRecord)
		}
		variants[record.Chr][record.Start] = append(
			variants[record.Chr][record.Start],
			record,
		)
	}

	if err := scanner.Err(); err != nil {
		return nil, "", err
	}

	return variants, sampleName, nil
}

func parseVariantRecord(fields []string) (*VariantRecord, error) {
	// Helper to parse int
	parseInt := func(s string) int {
		v, _ := strconv.Atoi(s)
		return v
	}

	// Helper to parse float
	parseFloat := func(s string) float64 {
		v, _ := strconv.ParseFloat(s, 64)
		return v
	}

	record := &VariantRecord{
		Sample:    fields[0],
		Gene:      fields[1],
		Chr:       fields[2],
		Start:     parseInt(fields[3]),
		End:       parseInt(fields[4]),
		Ref:       fields[5],
		Alt:       fields[6],
		Depth:     parseInt(fields[7]),
		AltDepth:  parseInt(fields[8]),
		RefFwd:    parseInt(fields[9]),
		RefRev:    parseInt(fields[10]),
		AltFwd:    parseInt(fields[11]),
		AltRev:    parseInt(fields[12]),
		Genotype:  fields[13],
		AF:        parseFloat(fields[14]),
		Bias:      fields[15],
		PMean:     parseFloat(fields[16]),
		PStd:      parseFloat(fields[17]),
		Qual:      parseFloat(fields[18]),
		QStd:      parseFloat(fields[19]),
		MapQ:      parseFloat(fields[20]),
		QRatio:    parseFloat(fields[21]),
		HiAF:      parseFloat(fields[22]),
		AdjAF:     parseFloat(fields[23]),
		NM:        parseFloat(fields[24]),
		MQ:        parseFloat(fields[25]),
		DupRate:   parseFloat(fields[26]),
		SV:        fields[27],
		Type:      fields[28],
		Status:    fields[29],
		PValue:    parseFloat(fields[30]),
		OddsRatio: parseFloat(fields[31]),
	}

	return record, nil
}

func processVariants(variants map[string]map[int][]*VariantRecord, config Config, writer io.Writer) {
	// Sort chromosomes
	chromosomes := sortChromosomes(variants)

	// Process each chromosome in order
	for _, chr := range chromosomes {
		processChromosome(chr, variants[chr], config, writer)
	}
}

func processChromosome(chr string, positions map[int][]*VariantRecord, config Config, writer io.Writer) {
	// Get sorted positions
	positionList := make([]int, 0, len(positions))
	for pos := range positions {
		positionList = append(positionList, pos)
	}
	sort.Ints(positionList)

	// Process each position
	for _, pos := range positionList {
		records := positions[pos]

		outputCount := 0
		for _, record := range records {
			// Apply filters
			filters := applyFilters(record, config)

			// Skip if marked for skipping
			if len(filters) == 1 && filters[0] == "SKIP" {
				continue
			}

			// Build and output VCF record
			vcfLine := buildVCFRecord(record, filters, config)

			// Handle -A flag (output all variants at same position)
			if config.AllVariants || outputCount == 0 {
				fmt.Fprint(writer, vcfLine)
				outputCount++
			}
		}
	}
}

func applyFilters(record *VariantRecord, config Config) []string {
	var filters []string

	// Calculate derived values
	refDepth := record.RefFwd + record.RefRev
	altDepth := record.AltFwd + record.AltRev
	totalDepth := refDepth + altDepth

	// Adjust odds ratio (handle Inf and invert if < 1)
	oddsRatio := record.OddsRatio
	if math.IsInf(oddsRatio, 1) || math.IsInf(oddsRatio, -1) {
		oddsRatio = 0
	} else if oddsRatio < 1 && oddsRatio > 0 {
		oddsRatio = 1 / oddsRatio
	}

	// Signal to Noise calculation
	sn := float64(altDepth) / (float64(totalDepth) + 0.5)

	// Check if this is an amplicon call (Bias == "1")
	isAmplicon := record.Bias == "1"

	// Filter 1: Minimum total depth
	if totalDepth < config.MinTotalDepth {
		if !(float64(record.AltDepth)*record.HiAF >= 0.5) {
			filters = append(filters, fmt.Sprintf("d%d", config.MinTotalDepth))
		}
	}

	// Filter 2: Minimum variant depth
	if record.AltDepth < config.MinVariantDepth {
		if !(float64(record.AltDepth)*record.HiAF >= 0.5) {
			filters = append(filters, fmt.Sprintf("v%d", config.MinVariantDepth))
		}
	}

	// Filter 3: Minimum allele frequency
	if record.AF < config.MinAF {
		filters = append(filters, fmt.Sprintf("f%s", formatFloat(config.MinAF)))
	}

	// Filter 4: Mean position in read
	if record.PMean < config.MinPMean {
		filters = append(filters, fmt.Sprintf("p%s", formatFloat(config.MinPMean)))
	}

	// Filter 5: Position standard deviation (if -P flag set)
	if config.FilterPStd && record.PStd == 0 && !isAmplicon && record.AF < 0.35 {
		filters = append(filters, "pSTD")
	}

	// Filter 6: Base quality
	if record.Qual < config.MinQual {
		filters = append(filters, fmt.Sprintf("q%s", formatFloat(config.MinQual)))
	}

	// Filter 7: Mapping quality
	if record.MapQ < config.MinMapQ && record.AF < 0.8 {
		filters = append(filters, fmt.Sprintf("Q%s", formatFloat(config.MinMapQ)))
	}

	// Filter 8: Signal to noise ratio
	if sn < config.MinSN {
		filters = append(filters, fmt.Sprintf("SN%s", formatFloat(config.MinSN)))
	}

	// Filter 9: Strand bias (Fisher's exact test)
	if record.PValue < 0.01 && oddsRatio >= 2.0 {
		filters = append(filters, "BIAS")
	}

	// Filter 10: Structural variant split read support
	if strings.Contains(record.Alt, "<") {
		splitReads := extractSplitReads(record.SV)
		if splitReads < config.MinSplitReads {
			return []string{"SKIP"}
		}
	}

	return filters
}

func buildVCFRecord(record *VariantRecord, filters []string, config Config) string {
	// Calculate values
	refDepth := record.RefFwd + record.RefRev
	altDepth := record.AltFwd + record.AltRev

	genotype := calculateGenotype(record, config)

	// Handle REF==ALT case (VCF spec requirement)
	ref := record.Ref
	alt := record.Alt
	if ref == alt {
		alt = "."
		genotype = "0/0"
	}

	// Build FILTER column
	filterStr := "PASS"
	if len(filters) > 0 {
		filterStr = strings.Join(filters, ";")
	}

	// Build INFO field
	info := buildINFO(record, config)

	// Build FORMAT field
	format := "GT:DP:VD:AD:AF:RD:ALD:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM:DUPRATE"

	// Build SAMPLE field
	sample := buildSampleData(record, genotype, refDepth, altDepth)

	// Construct VCF line
	return fmt.Sprintf("%s\t%d\t.\t%s\t%s\t.\t%s\t%s\t%s\t%s\n",
		record.Chr,
		record.Start,
		ref,
		alt,
		filterStr,
		info,
		format,
		sample,
	)
}

func calculateGenotype(record *VariantRecord, config Config) string {
	if record.Ref == record.Alt {
		return "0/0"
	}

	af := record.AF

	if (1 - af) < config.GTFreq {
		return "1/1" // Homozygous alt
	} else if af >= 0.5 {
		return "1/0" // Heterozygous (alt/ref)
	} else if af >= config.MinAF {
		return "0/1" // Heterozygous (ref/alt)
	}

	return "0/0" // Reference
}

func buildINFO(record *VariantRecord, config Config) string {
	var parts []string

	sampleNoSpace := strings.ReplaceAll(record.Sample, " ", "_")
	parts = append(parts, fmt.Sprintf("SAMPLE=%s", sampleNoSpace))
	parts = append(parts, fmt.Sprintf("TYPE=%s", record.Type))
	parts = append(parts, fmt.Sprintf("DP=%d", record.Depth))

	if config.PrintEndTag {
		parts = append(parts, fmt.Sprintf("END=%d", record.End))
	}

	parts = append(parts, fmt.Sprintf("VD=%d", record.AltDepth))
	parts = append(parts, fmt.Sprintf("AF=%.4f", record.AF))
	parts = append(parts, fmt.Sprintf("BIAS=%s", record.Bias))
	parts = append(parts, fmt.Sprintf("PMEAN=%.1f", record.PMean))
	parts = append(parts, fmt.Sprintf("PSTD=%.1f", record.PStd))
	parts = append(parts, fmt.Sprintf("QUAL=%.1f", record.Qual))
	parts = append(parts, fmt.Sprintf("QSTD=%.1f", record.QStd))
	parts = append(parts, fmt.Sprintf("SBF=%d:%d:%d:%d",
		record.RefFwd, record.RefRev, record.AltFwd, record.AltRev))

	// Adjust odds ratio for INFO
	oddsRatio := record.OddsRatio
	if math.IsInf(oddsRatio, 1) || math.IsInf(oddsRatio, -1) {
		oddsRatio = 0
	} else if oddsRatio < 1 && oddsRatio > 0 {
		oddsRatio = 1 / oddsRatio
	}
	parts = append(parts, fmt.Sprintf("ODDRATIO=%.4f", oddsRatio))

	parts = append(parts, fmt.Sprintf("MQ=%.1f", record.MapQ))
	parts = append(parts, fmt.Sprintf("SN=%.4f", record.QRatio))
	parts = append(parts, fmt.Sprintf("HIAF=%.4f", record.HiAF))
	parts = append(parts, fmt.Sprintf("ADJAF=%.4f", record.AdjAF))
	parts = append(parts, fmt.Sprintf("NM=%.1f", record.NM))
	parts = append(parts, fmt.Sprintf("DUPRATE=%.4f", record.DupRate))

	// Add structural variant info if present
	if strings.Contains(record.Alt, "<") {
		svInfo := buildSVInfo(record)
		if svInfo != "" {
			parts = append(parts, svInfo)
		}
	}

	return strings.Join(parts, ";")
}

func buildSVInfo(record *VariantRecord) string {
	// Extract SV type from ALT field
	svType := strings.Trim(record.Alt, "<>")
	svLen := record.End - record.Start
	if strings.Contains(svType, "INV") {
		svLen++
	}

	parts := []string{
		fmt.Sprintf("SVTYPE=%s", svType),
		fmt.Sprintf("SVLEN=%d", svLen),
	}

	// Add split read and spanning pair info if available
	if record.SV != "" {
		if splitReads := extractSplitReads(record.SV); splitReads > 0 {
			spanPairs := extractSpanPairs(record.SV)
			parts = append(parts,
				fmt.Sprintf("SPLITREAD=%d", splitReads),
				fmt.Sprintf("SPANPAIR=%d", spanPairs),
			)
		}
	}

	return strings.Join(parts, ";")
}

func buildSampleData(record *VariantRecord, gt string, refDepth, altDepth int) string {
	// Adjust odds ratio
	oddsRatio := record.OddsRatio
	if math.IsInf(oddsRatio, 1) || math.IsInf(oddsRatio, -1) {
		oddsRatio = 0
	} else if oddsRatio < 1 && oddsRatio > 0 {
		oddsRatio = 1 / oddsRatio
	}

	sbf := fmt.Sprintf("%d:%d:%d:%d",
		record.RefFwd, record.RefRev,
		record.AltFwd, record.AltRev)

	return fmt.Sprintf("%s:%d:%d:%d,%d:%.4f:%d:%d:%s:%.1f:%.1f:%.1f:%.1f:%s:%.4f:%.1f:%.4f:%.4f:%.4f:%.1f:%.4f",
		gt,
		record.Depth,
		record.AltDepth,
		refDepth, altDepth,
		record.AF,
		refDepth,
		altDepth,
		record.Bias,
		record.PMean,
		record.PStd,
		record.Qual,
		record.QStd,
		sbf,
		oddsRatio,
		record.MapQ,
		record.QRatio,
		record.HiAF,
		record.AdjAF,
		record.NM,
		record.DupRate,
	)
}

func extractSplitReads(svField string) int {
	if svField == "" {
		return 0
	}
	parts := strings.Split(svField, "-")
	if len(parts) > 0 {
		sr, _ := strconv.Atoi(parts[0])
		return sr
	}
	return 0
}

func extractSpanPairs(svField string) int {
	if svField == "" {
		return 0
	}
	parts := strings.Split(svField, "-")
	if len(parts) > 1 {
		sp, _ := strconv.Atoi(parts[1])
		return sp
	}
	return 0
}

func formatFloat(f float64) string {
	s := fmt.Sprintf("%.4f", f)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	return s
}

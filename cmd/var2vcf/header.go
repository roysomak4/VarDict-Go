package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

func printVCFHeader(config Config, sampleName string, writer io.Writer) {
	// File format
	fmt.Fprintf(writer, "##fileformat=VCFv4.2\n")
	fmt.Fprintf(writer, "##source=VarDict_Go_v%s\n", Version)

	// Reference
	if config.RefPath != "" {
		fmt.Fprintf(writer, "##reference=%s\n", config.RefPath)
	}

	// Contigs from BED file
	printContigs(config.BedPath, writer)

	// INFO field definitions
	printINFOHeaders(writer)

	// FORMAT field definitions
	printFORMATHeaders(writer)

	// FILTER definitions
	printFILTERHeaders(config, writer)

	// Column header line
	sampleNoSpace := strings.ReplaceAll(sampleName, " ", "_")
	fmt.Fprintf(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sampleNoSpace)
}

func printContigs(bedPath string, writer io.Writer) {
	if bedPath == "" {
		return
	}

	file, err := os.Open(bedPath)
	if err != nil {
		return
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) >= 3 {
			chr := fields[0]
			end := fields[2]
			fmt.Fprintf(writer, "##contig=<ID=%s,length=%s>\n", chr, end)
		}
	}
}

func printINFOHeaders(writer io.Writer) {
	infoFields := []struct {
		id          string
		number      string
		typ         string
		description string
	}{
		{"SAMPLE", "1", "String", "Sample name (with whitespace translated to underscores)"},
		{"TYPE", "1", "String", "Variant Type: SNV Insertion Deletion Complex"},
		{"DP", "1", "Integer", "Total Depth"},
		{"END", "1", "Integer", "Chr End Position"},
		{"VD", "1", "Integer", "Variant Depth"},
		{"AF", "A", "Float", "Allele Frequency"},
		{"BIAS", "1", "String", "Strand Bias Info"},
		{"PMEAN", "1", "Float", "Mean position in reads"},
		{"PSTD", "1", "Float", "Position STD in reads"},
		{"QUAL", "1", "Float", "Mean base quality"},
		{"QSTD", "1", "Float", "Base quality STD"},
		{"SBF", "1", "String", "Strand Bias Fisher: RefFor:RefRev:AltFor:AltRev"},
		{"ODDRATIO", "1", "Float", "Strand Bias Odds Ratio"},
		{"MQ", "1", "Float", "Mean mapping quality"},
		{"SN", "1", "Float", "Signal to noise ratio"},
		{"HIAF", "1", "Float", "High quality allele frequency"},
		{"ADJAF", "1", "Float", "Adjusted allele frequency"},
		{"NM", "1", "Float", "Mean mismatches in reads"},
		{"DUPRATE", "1", "Float", "Duplication rate"},
		{"SVTYPE", "1", "String", "Structural variant type"},
		{"SVLEN", "1", "Integer", "Structural variant length"},
		{"SPLITREAD", "1", "Integer", "Number of split reads supporting SV"},
		{"SPANPAIR", "1", "Integer", "Number of spanning pairs supporting SV"},
	}

	for _, field := range infoFields {
		fmt.Fprintf(writer, "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
			field.id, field.number, field.typ, field.description)
	}
}

func printFORMATHeaders(writer io.Writer) {
	formatFields := []struct {
		id          string
		number      string
		typ         string
		description string
	}{
		{"GT", "1", "String", "Genotype"},
		{"DP", "1", "Integer", "Total Depth"},
		{"VD", "1", "Integer", "Variant Depth"},
		{"AD", "R", "Integer", "Allelic depths for the ref and alt alleles"},
		{"AF", "A", "Float", "Allele Frequency"},
		{"RD", "1", "Integer", "Reference Depth"},
		{"ALD", "1", "Integer", "Alternate allele depth"},
		{"BIAS", "1", "String", "Strand Bias"},
		{"PMEAN", "1", "Float", "Mean position in reads"},
		{"PSTD", "1", "Float", "Position STD in reads"},
		{"QUAL", "1", "Float", "Mean base quality"},
		{"QSTD", "1", "Float", "Base quality STD"},
		{"SBF", "1", "String", "Strand Bias Fisher"},
		{"ODDRATIO", "1", "Float", "Strand Bias Odds Ratio"},
		{"MQ", "1", "Float", "Mean mapping quality"},
		{"SN", "1", "Float", "Signal to noise"},
		{"HIAF", "1", "Float", "High quality allele frequency"},
		{"ADJAF", "1", "Float", "Adjusted allele frequency"},
		{"NM", "1", "Float", "Mean mismatches"},
		{"DUPRATE", "1", "Float", "Duplication rate"},
	}

	for _, field := range formatFields {
		fmt.Fprintf(writer, "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
			field.id, field.number, field.typ, field.description)
	}
}

func printFILTERHeaders(config Config, writer io.Writer) {
	filters := []struct {
		id          string
		description string
	}{
		{"PASS", "Passed all filters"},
		{fmt.Sprintf("d%d", config.MinTotalDepth), fmt.Sprintf("Total depth < %d", config.MinTotalDepth)},
		{fmt.Sprintf("v%d", config.MinVariantDepth), fmt.Sprintf("Variant depth < %d", config.MinVariantDepth)},
		{fmt.Sprintf("f%s", formatFloat(config.MinAF)), fmt.Sprintf("Allele frequency < %s", formatFloat(config.MinAF))},
		{fmt.Sprintf("p%s", formatFloat(config.MinPMean)), fmt.Sprintf("Mean position < %s", formatFloat(config.MinPMean))},
		{"pSTD", "Position standard deviation = 0"},
		{fmt.Sprintf("q%s", formatFloat(config.MinQual)), fmt.Sprintf("Mean base quality < %s", formatFloat(config.MinQual))},
		{fmt.Sprintf("Q%s", formatFloat(config.MinMapQ)), fmt.Sprintf("Mean mapping quality < %s", formatFloat(config.MinMapQ))},
		{fmt.Sprintf("SN%s", formatFloat(config.MinSN)), fmt.Sprintf("Signal to noise < %s", formatFloat(config.MinSN))},
		{"BIAS", "Strand bias detected (p-value < 0.01 and odds ratio >= 2)"},
	}

	for _, filter := range filters {
		fmt.Fprintf(writer, "##FILTER=<ID=%s,Description=\"%s\">\n",
			filter.id, filter.description)
	}
}

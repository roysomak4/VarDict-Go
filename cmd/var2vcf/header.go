package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

func printVCFHeader(config Config, sampleName string, writer io.Writer) {
	fmt.Fprintf(writer, "##fileformat=VCFv4.2\n")
	fmt.Fprintf(writer, "##source=VarDict_Go_v%s\n", Version)

	if config.RefPath != "" {
		fmt.Fprintf(writer, "##reference=%s\n", config.RefPath)
	}

	printContigs(config.BedPath, writer)
	printINFOHeaders(writer)
	printFILTERHeaders(config, writer)
	printFORMATHeaders(writer)

	sampleNoSpace := strings.ReplaceAll(sampleName, " ", "_")
	fmt.Fprintf(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sampleNoSpace)
}

func printContigs(bedPath string, writer io.Writer) {
	if bedPath == "" {
		return
	}
	f, err := os.Open(bedPath)
	if err != nil {
		return
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) >= 3 {
			fmt.Fprintf(writer, "##contig=<ID=%s,length=%s>\n", fields[0], fields[2])
		}
	}
}

func printINFOHeaders(writer io.Writer) {
	fields := []struct{ id, number, typ, desc string }{
		{"SAMPLE", "1", "String", "Sample name (with whitespace translated to underscores)"},
		{"TYPE", "1", "String", "Variant Type: SNV Insertion Deletion Complex"},
		{"DP", "1", "Integer", "Total Depth"},
		{"END", "1", "Integer", "Chr End Position"},
		{"VD", "1", "Integer", "Variant Depth"},
		{"AF", "A", "Float", "Allele Frequency"},
		{"BIAS", "1", "String", "Strand Bias Info"},
		{"REFBIAS", "1", "String", "Reference depth by strand"},
		{"VARBIAS", "1", "String", "Variant depth by strand"},
		{"PMEAN", "1", "Float", "The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call"},
		{"PSTD", "1", "Float", "Position STD in reads"},
		{"QUAL", "1", "Float", "Mean quality score in reads"},
		{"QSTD", "1", "Float", "Quality score STD in reads"},
		{"SBF", "1", "Float", "Strand Bias Fisher p-value"},
		{"ODDRATIO", "1", "Float", "Strand Bias Odds ratio"},
		{"MQ", "1", "Float", "Mean Mapping Quality"},
		{"SN", "1", "Float", "Signal to noise"},
		{"HIAF", "1", "Float", "Allele frequency using only high quality bases"},
		{"ADJAF", "1", "Float", "Adjusted AF for indels due to local realignment"},
		{"SHIFT3", "1", "Integer", "No. of bases to be shifted to 3 prime for deletions due to alternative alignment"},
		{"MSI", "1", "Float", "MicroSatellite. > 1 indicates MSI"},
		{"MSILEN", "1", "Float", "MicroSatellite unit length in bp"},
		{"NM", "1", "Float", "Mean mismatches in reads"},
		{"LSEQ", "1", "String", "5' flanking seq"},
		{"RSEQ", "1", "String", "3' flanking seq"},
		{"GDAMP", "1", "Integer", "No. of amplicons supporting variant"},
		{"TLAMP", "1", "Integer", "Total of amplicons covering variant"},
		{"NCAMP", "1", "Integer", "No. of amplicons don't work"},
		{"AMPFLAG", "1", "Integer", "Top variant in amplicons don't match"},
		{"HICNT", "1", "Integer", "High quality variant reads"},
		{"HICOV", "1", "Integer", "High quality total reads"},
		{"SPLITREAD", "1", "Integer", "No. of split reads supporting SV"},
		{"SPANPAIR", "1", "Integer", "No. of pairs supporting SV"},
		{"SVTYPE", "1", "String", "SV type: INV DUP DEL INS FUS"},
		{"SVLEN", "1", "Integer", "The length of SV in bp"},
		{"DUPRATE", "1", "Float", "Duplication rate in fraction"},
	}
	for _, f := range fields {
		fmt.Fprintf(writer, "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
			f.id, f.number, f.typ, f.desc)
	}
}

func printFILTERHeaders(config Config, writer io.Writer) {
	filters := []struct{ id, desc string }{
		{fmt.Sprintf("q%g", config.MinQual), fmt.Sprintf("Mean Base Quality Below %g", config.MinQual)},
		{fmt.Sprintf("Q%g", config.MinMapQ), fmt.Sprintf("Mean Mapping Quality Below %g", config.MinMapQ)},
		{fmt.Sprintf("p%g", config.MinPMean), fmt.Sprintf("Mean Position in Reads Less than %g", config.MinPMean)},
		{fmt.Sprintf("SN%g", config.MinSN), fmt.Sprintf("Signal to Noise Less than %g", config.MinSN)},
		{"Bias", "Strand Bias"},
		{"pSTD", "Position in Reads has STD of 0"},
		{fmt.Sprintf("d%d", config.MinTotalDepth), fmt.Sprintf("Total Depth < %d", config.MinTotalDepth)},
		{fmt.Sprintf("v%d", config.MinVarDepth), fmt.Sprintf("Var Depth < %d", config.MinVarDepth)},
		{fmt.Sprintf("f%g", config.MinAF), fmt.Sprintf("Allele frequency < %g", config.MinAF)},
		{fmt.Sprintf("MSI%d", config.MaxMSI), fmt.Sprintf("Variant in MSI region with %d non-monomer MSI or 13 monomer MSI", config.MaxMSI)},
		{fmt.Sprintf("NM%g", config.MaxNM), fmt.Sprintf("Mean mismatches in reads >= %g, thus likely false positive", config.MaxNM)},
		{"InGap", "The variant is in the deletion gap, thus likely false positive"},
		{"InIns", "The variant is adjacent to an insertion variant"},
		{fmt.Sprintf("Cluster%dbp", config.ClusterBP), fmt.Sprintf("Two variants are within %d bp", config.ClusterBP)},
		{"LongMSI", "The somatic variant is flanked by long A/T (>=14)"},
		{"AMPBIAS", "Indicate the variant has amplicon bias."},
	}
	for _, f := range filters {
		fmt.Fprintf(writer, "##FILTER=<ID=%s,Description=\"%s\">\n", f.id, f.desc)
	}
}

func printFORMATHeaders(writer io.Writer) {
	fields := []struct{ id, number, typ, desc string }{
		{"GT", "1", "String", "Genotype"},
		{"DP", "1", "Integer", "Total Depth"},
		{"VD", "1", "Integer", "Variant Depth"},
		{"AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed"},
		{"AF", "A", "Float", "Allele Frequency"},
		{"RD", "2", "Integer", "Reference forward, reverse reads"},
		{"ALD", "2", "Integer", "Variant forward, reverse reads"},
	}
	for _, f := range fields {
		fmt.Fprintf(writer, "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
			f.id, f.number, f.typ, f.desc)
	}
}

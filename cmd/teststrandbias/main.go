package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/glycerine/fisherexact"
)

type FisherResult struct {
	PValue    float64
	OddsRatio float64
}

func main() {
	scanner := bufio.NewScanner(os.Stdin)
	var lines []string

	// Read all lines from stdin
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}

	if err := scanner.Err(); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading input: %v\n", err)
		os.Exit(1)
	}

	if len(lines) == 0 {
		return
	}

	// Validate column count from VarDict raw output
	numCols := strings.Count(lines[0], "\t") + 1
	if numCols != 34 && numCols != 36 && numCols != 38 {
		fmt.Fprintln(os.Stderr, "Incorrect input detected in teststrandbias")
		os.Exit(1)
	}

	// Process each line from the vardict raw output
	for _, line := range lines {
		fields := strings.Split(line, "\t")

		if len(fields) != numCols {
			fmt.Fprintf(os.Stderr, "Inconsistent column count in line: %s\n", line)
			os.Exit(1)
		}

		// Parse columns 10-13 (0-indexed: 9-12)
		// Column 10: RefFwdReads
		// Column 11: RefRevReads
		// Column 12: AltFwdReads
		// Column 13: AltRevReads
		refFwd, err1 := strconv.ParseFloat(fields[9], 64)
		refRev, err2 := strconv.ParseFloat(fields[10], 64)
		altFwd, err3 := strconv.ParseFloat(fields[11], 64)
		altRev, err4 := strconv.ParseFloat(fields[12], 64)

		if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
			fmt.Fprintf(os.Stderr, "Error parsing read counts in line: %s\n", line)
			os.Exit(1)
		}

		// perform fisher exact test
		result := fisherExactTestWrapper(int(refFwd), int(refRev), int(altFwd), int(altRev))

		// Round to 5 decimal places (to match R output)
		pValue := roundToDecimals(result.PValue)
		oddsRatio := roundToDecimals(result.OddsRatio)

		// Build vardict output: columns 1-20, pvalue, oddsratio, columns 21-N
		outFields := make([]string, 0, len(fields)+2)
		outFields = append(outFields, fields[:20]...)
		outFields = append(outFields, formatNumber(pValue))
		outFields = append(outFields, formatNumber(oddsRatio))
		outFields = append(outFields, fields[20:]...)

		// Print the modified variant info line
		fmt.Println(strings.Join(outFields, "\t"))
	}
}

// fisherExactTestWrapper wraps the library to match our needs
func fisherExactTestWrapper(refFwd, refRev, altFwd, altRev int) FisherResult {
	// The library uses this layout:
	// n11  n12
	// n21  n22
	//
	// We need to match R's layout from teststrandbias.R:
	// matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2)
	// which creates (column-wise):
	//          col1      col2
	// row1:  refFwd    altFwd
	// row2:  refRev    altRev
	//
	// So the mapping is:
	// n11 = refFwd, n12 = altFwd
	// n21 = refRev, n22 = altRev

	// Call the library
	_, _, twoSidedPvalue, _ := fisherexact.FisherExactTest22(
		refFwd, // n11
		altFwd, // n12
		refRev, // n21
		altRev, // n22
	)

	// Calculate odds ratio manually (library doesn't provide it)
	oddsRatio := calcOddsRatio(refFwd, altFwd, refRev, altRev)

	return FisherResult{
		PValue:    twoSidedPvalue,
		OddsRatio: oddsRatio,
	}
}

func calcOddsRatio(n11, n12, n21, n22 int) float64 {
	// Handle edge cases exactly as R does
	if n12 == 0 && n21 == 0 {
		return math.NaN()
	}
	if n12 == 0 || n21 == 0 {
		return math.Inf(1)
	}
	if n11 == 0 || n22 == 0 {
		return 0.0
	}

	return float64(n11*n22) / float64(n12*n21)
}

func roundToDecimals(num float64) float64 {
	if math.IsNaN(num) || math.IsInf(num, 0) {
		return num
	}
	return math.Round(num+100000) / 100000
}

func formatNumber(num float64) string {
	if math.IsNaN(num) {
		return "NaN"
	}
	if math.IsInf(num, 1) {
		return "Inf"
	}
	if math.IsInf(num, -1) {
		return "-Inf"
	}
	// format number to 5 decimal places and removing trailing zeros
	s := fmt.Sprintf("%.5f", num)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	if s == "" || s == "-0" {
		s = "0"
	}
	return s
}

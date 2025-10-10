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

// FisherResult holds the results for our VarDict use case
type FisherResult struct {
	PValue    float64
	OddsRatio float64
}

func main() {
	scanner := bufio.NewScanner(os.Stdin)
	var lines []string

	// Read all input from stdin
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

	// Validate column count
	numCols := strings.Count(lines[0], "\t") + 1
	if numCols != 34 && numCols != 36 && numCols != 38 {
		fmt.Fprintln(os.Stderr, "Incorrect input detected in teststrandbias")
		os.Exit(1)
	}

	// Process each line
	for _, line := range lines {
		fields := strings.Split(line, "\t")

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
			fmt.Fprintf(os.Stderr, "Error parsing strand counts\n")
			os.Exit(1)
		}

		// Perform Fisher's exact test using library
		result := fisherExactTestWrapper(int(refFwd), int(refRev), int(altFwd), int(altRev))

		// Round to 5 decimal places (matching R)
		pvalue := roundTo5Decimals(result.PValue)
		oddsRatio := roundTo5Decimals(result.OddsRatio)

		// Build output: columns 1-20, pvalue, oddsratio, columns 21-N
		outputFields := make([]string, 0, len(fields)+2)
		outputFields = append(outputFields, fields[0:20]...)
		outputFields = append(outputFields, formatNumber(pvalue))
		outputFields = append(outputFields, formatNumber(oddsRatio))
		if len(fields) > 20 {
			outputFields = append(outputFields, fields[20:]...)
		}

		fmt.Println(strings.Join(outputFields, "\t"))
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
	oddsRatio := calculateOddsRatio(refFwd, altFwd, refRev, altRev)

	return FisherResult{
		PValue:    twoSidedPvalue,
		OddsRatio: oddsRatio,
	}
}

// calculateOddsRatio computes the odds ratio matching R's fisher.test()
// CRITICAL: Order of checks must match R's exact behavior!
func calculateOddsRatio(n11, n12, n21, n22 int) float64 {
	// R's fisher.test checks numerator FIRST, then denominator
	// OR = (n11 * n22) / (n12 * n21)

	// If numerator is zero, OR = 0 (even if denominator is also zero)
	if n11 == 0 || n22 == 0 {
		return 0.0
	}

	// If denominator is zero (but numerator isn't), OR = Inf
	if n12 == 0 || n21 == 0 {
		return math.Inf(1)
	}

	// Standard odds ratio calculation
	return float64(n11*n22) / float64(n12*n21)
}

// roundTo5Decimals rounds a float to 5 decimal places
func roundTo5Decimals(x float64) float64 {
	return math.Round(x*100000) / 100000
}

// formatNumber formats a number without scientific notation
func formatNumber(x float64) string {
	if math.IsNaN(x) {
		return "NaN"
	}
	if math.IsInf(x, 1) {
		return "Inf"
	}
	if math.IsInf(x, -1) {
		return "-Inf"
	}

	// Format with up to 5 decimal places, removing trailing zeros
	s := fmt.Sprintf("%.5f", x)

	// Remove trailing zeros after decimal point
	if strings.Contains(s, ".") {
		s = strings.TrimRight(s, "0")
		s = strings.TrimRight(s, ".")
	}

	return s
}

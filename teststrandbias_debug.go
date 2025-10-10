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
	fmt.Fprintf(os.Stderr, "DEBUG: Number of columns: %d\n", numCols)

	if numCols != 34 && numCols != 36 && numCols != 38 {
		fmt.Fprintln(os.Stderr, "Incorrect input detected in teststrandbias")
		os.Exit(1)
	}

	// Process only first 3 lines for debugging
	maxLines := len(lines)
	if maxLines > 3 {
		maxLines = 3
	}

	for lineNum, line := range lines[:maxLines] {
		fields := strings.Split(line, "\t")

		fmt.Fprintf(os.Stderr, "\n=== Line %d ===\n", lineNum+1)
		fmt.Fprintf(os.Stderr, "Total fields: %d\n", len(fields))

		// Print first 20 columns to see structure
		fmt.Fprintf(os.Stderr, "First 15 fields:\n")
		for i := 0; i < 15 && i < len(fields); i++ {
			fmt.Fprintf(os.Stderr, "  [%d]: %s\n", i, fields[i])
		}

		// Parse columns 9-12 (0-indexed) which should be columns 10-13 in VarDict (1-indexed)
		fmt.Fprintf(os.Stderr, "\nStrand count columns (0-indexed 9-12):\n")
		if len(fields) > 12 {
			fmt.Fprintf(os.Stderr, "  fields[9] (RefFwd): %s\n", fields[9])
			fmt.Fprintf(os.Stderr, "  fields[10] (RefRev): %s\n", fields[10])
			fmt.Fprintf(os.Stderr, "  fields[11] (AltFwd): %s\n", fields[11])
			fmt.Fprintf(os.Stderr, "  fields[12] (AltRev): %s\n", fields[12])
		}

		// Try parsing
		refFwd, err1 := strconv.ParseFloat(fields[9], 64)
		refRev, err2 := strconv.ParseFloat(fields[10], 64)
		altFwd, err3 := strconv.ParseFloat(fields[11], 64)
		altRev, err4 := strconv.ParseFloat(fields[12], 64)

		if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
			fmt.Fprintf(os.Stderr, "ERROR parsing: %v %v %v %v\n", err1, err2, err3, err4)
			continue
		}

		fmt.Fprintf(os.Stderr, "\nParsed values:\n")
		fmt.Fprintf(os.Stderr, "  RefFwd: %.0f\n", refFwd)
		fmt.Fprintf(os.Stderr, "  RefRev: %.0f\n", refRev)
		fmt.Fprintf(os.Stderr, "  AltFwd: %.0f\n", altFwd)
		fmt.Fprintf(os.Stderr, "  AltRev: %.0f\n", altRev)

		// Calculate Fisher's test
		_, _, pvalue, _ := fisherexact.FisherExactTest22(
			int(refFwd),
			int(altFwd),
			int(refRev),
			int(altRev),
		)

		// Calculate odds ratio
		or := calculateOddsRatio(int(refFwd), int(altFwd), int(refRev), int(altRev))

		fmt.Fprintf(os.Stderr, "\nResults:\n")
		fmt.Fprintf(os.Stderr, "  P-value: %.5f\n", pvalue)
		fmt.Fprintf(os.Stderr, "  Odds Ratio: %.5f\n", or)

		// Also show what R would get with same values
		fmt.Fprintf(os.Stderr, "\nR command would be:\n")
		fmt.Fprintf(os.Stderr, "  fisher.test(matrix(c(%.0f, %.0f, %.0f, %.0f), nrow=2))\n",
			refFwd, refRev, altFwd, altRev)
	}
}

func calculateOddsRatio(n11, n12, n21, n22 int) float64 {
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

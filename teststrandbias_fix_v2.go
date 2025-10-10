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

	numCols := strings.Count(lines[0], "\t") + 1
	if numCols != 34 && numCols != 36 && numCols != 38 {
		fmt.Fprintln(os.Stderr, "Incorrect input detected in teststrandbias")
		os.Exit(1)
	}

	for _, line := range lines {
		fields := strings.Split(line, "\t")

		refFwd, err1 := strconv.ParseFloat(fields[9], 64)
		refRev, err2 := strconv.ParseFloat(fields[10], 64)
		altFwd, err3 := strconv.ParseFloat(fields[11], 64)
		altRev, err4 := strconv.ParseFloat(fields[12], 64)

		if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
			fmt.Fprintf(os.Stderr, "Error parsing strand counts\n")
			os.Exit(1)
		}

		result := fisherExactTestWrapper(int(refFwd), int(refRev), int(altFwd), int(altRev))

		pvalue := roundTo5Decimals(result.PValue)
		oddsRatio := roundTo5Decimals(result.OddsRatio)

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

func fisherExactTestWrapper(refFwd, refRev, altFwd, altRev int) FisherResult {
	_, _, twoSidedPvalue, _ := fisherexact.FisherExactTest22(
		refFwd, altFwd, refRev, altRev,
	)

	// Calculate conditional MLE of odds ratio (matches R)
	oddsRatio := calculateConditionalMLEOddsRatio(refFwd, altFwd, refRev, altRev)

	return FisherResult{
		PValue:    twoSidedPvalue,
		OddsRatio: oddsRatio,
	}
}

// calculateConditionalMLEOddsRatio computes the conditional MLE matching R's fisher.test()
func calculateConditionalMLEOddsRatio(n11, n12, n21, n22 int) float64 {
	// Handle edge cases first
	if n11 == 0 || n22 == 0 {
		return 0.0
	}
	if n12 == 0 || n21 == 0 {
		return math.Inf(1)
	}

	// For conditional MLE, we need to find the odds ratio that maximizes
	// the probability of observing this table given the marginals
	// This is done via Newton-Raphson or similar numerical method

	// Calculate marginals
	row1 := n11 + n12
	row2 := n21 + n22
	col1 := n11 + n21

	// Start with simple odds ratio as initial guess
	or := float64(n11*n22) / float64(n12*n21)

	// Newton-Raphson iteration to find conditional MLE
	for iter := 0; iter < 50; iter++ {
		// Calculate expected value of n11 given current OR
		expN11 := expectedN11(or, row1, row2, col1)

		// Calculate derivative
		deriv := derivativeExpN11(or, row1, row2, col1)

		if math.Abs(deriv) < 1e-10 {
			break
		}

		// Newton-Raphson update
		delta := (float64(n11) - expN11) / deriv
		orNew := or + delta

		// Ensure OR stays positive
		if orNew <= 0 {
			orNew = or * 0.5
		}

		// Check convergence
		if math.Abs(orNew-or) < 1e-8*or {
			break
		}

		or = orNew
	}

	return or
}

// expectedN11 calculates E[n11] for given odds ratio under hypergeometric distribution
func expectedN11(or float64, row1, row2, col1 int) float64 {
	n := row1 + row2
	minVal := max(0, col1-row2)
	maxVal := min(row1, col1)

	var numerator, denominator float64

	for k := minVal; k <= maxVal; k++ {
		// Calculate hypergeometric probability with given OR
		logProb := logHypergeometricOR(k, row1, row2, col1, n, or)
		prob := math.Exp(logProb)

		numerator += float64(k) * prob
		denominator += prob
	}

	if denominator == 0 {
		return float64(col1*row1) / float64(n)
	}

	return numerator / denominator
}

// derivativeExpN11 calculates d/d(OR) E[n11]
func derivativeExpN11(or float64, row1, row2, col1 int) float64 {
	n := row1 + row2
	minVal := max(0, col1-row2)
	maxVal := min(row1, col1)

	var numerator, denominator float64
	expN11 := expectedN11(or, row1, row2, col1)

	for k := minVal; k <= maxVal; k++ {
		logProb := logHypergeometricOR(k, row1, row2, col1, n, or)
		prob := math.Exp(logProb)

		// Derivative involves (k - E[n11])^2
		diff := float64(k) - expN11
		numerator += diff * diff * prob
		denominator += prob
	}

	if denominator == 0 || or == 0 {
		return 0
	}

	return numerator / (denominator * or)
}

// logHypergeometricOR calculates log probability for hypergeometric with odds ratio
func logHypergeometricOR(k, row1, row2, col1, n int, or float64) float64 {
	// Standard hypergeometric probability
	logProb := logChoose(row1, k) + logChoose(row2, col1-k) - logChoose(n, col1)

	// Adjust by odds ratio
	if or != 1.0 && k > 0 {
		logProb += float64(k) * math.Log(or)
	}

	return logProb
}

// logChoose calculates log(n choose k)
func logChoose(n, k int) float64 {
	if k < 0 || k > n {
		return math.Inf(-1)
	}
	if k == 0 || k == n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	return logFactorial(n) - logFactorial(k) - logFactorial(n-k)
}

// logFactorial calculates log(n!)
func logFactorial(n int) float64 {
	if n < 0 {
		return math.NaN()
	}
	if n <= 1 {
		return 0
	}

	// Use lgamma for accuracy
	return lgamma(float64(n + 1))
}

// lgamma approximates log(Gamma(x)) using Lanczos approximation
func lgamma(x float64) float64 {
	if x < 0 {
		return math.NaN()
	}
	if x < 0.5 {
		return math.Log(math.Pi) - math.Log(math.Sin(math.Pi*x)) - lgamma(1-x)
	}

	// Lanczos approximation
	coef := []float64{
		0.99999999999980993, 676.5203681218851, -1259.1392167224028,
		771.32342877765313, -176.61502916214059, 12.507343278686905,
		-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
	}

	z := x
	z -= 1
	base := z + 7.5
	sum := coef[0]
	for i := 1; i < len(coef); i++ {
		sum += coef[i] / (z + float64(i))
	}

	return math.Log(math.Sqrt(2*math.Pi)) + math.Log(sum) - base + math.Log(base)*(z+0.5)
}

func roundTo5Decimals(x float64) float64 {
	return math.Round(x*100000) / 100000
}

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

	s := fmt.Sprintf("%.5f", x)
	if strings.Contains(s, ".") {
		s = strings.TrimRight(s, "0")
		s = strings.TrimRight(s, ".")
	}
	return s
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

package main

import (
	"sort"
	"strconv"
	"strings"
)

// sortChromosomes returns chromosomes in the correct order:
// 1. Numeric chromosomes (1, 2, ..., 22 or chr1, chr2, ..., chr22)
// 2. Sex chromosomes (X, Y or chrX, chrY)
// 3. Mitochondrial (MT, M, chrM, chrMT)
// 4. All other chromosomes (alphabetically)
func sortChromosomes(variants map[string]map[int][]*VariantRecord) []string {
	// Separate chromosomes into categories
	var numeric []int
	var other []string
	hasPrefix := false

	// First pass: determine if we have "chr" prefix
	for chr := range variants {
		if strings.HasPrefix(chr, "chr") {
			hasPrefix = true
			break
		}
	}

	// Second pass: categorize chromosomes
	for chr := range variants {
		chrNum := strings.TrimPrefix(chr, "chr")

		// Try to parse as integer
		if num, err := strconv.Atoi(chrNum); err == nil {
			numeric = append(numeric, num)
		} else {
			// Not a numeric chromosome
			other = append(other, chr)
		}
	}

	// Sort numeric chromosomes
	sort.Ints(numeric)

	// Build result with numeric chromosomes first
	result := make([]string, 0, len(variants))
	for _, num := range numeric {
		if hasPrefix {
			result = append(result, "chr"+strconv.Itoa(num))
		} else {
			result = append(result, strconv.Itoa(num))
		}
	}

	// Add sex chromosomes in order: X, Y
	sexChroms := []string{"X", "Y"}
	for _, sex := range sexChroms {
		testChr := sex
		if hasPrefix {
			testChr = "chr" + sex
		}
		if _, exists := variants[testChr]; exists {
			result = append(result, testChr)
			// Remove from other list if present
			other = removeString(other, testChr)
		}
	}

	// Add mitochondrial chromosomes: MT, M
	mtChroms := []string{"MT", "M"}
	for _, mt := range mtChroms {
		testChr := mt
		if hasPrefix {
			testChr = "chr" + mt
		}
		if _, exists := variants[testChr]; exists {
			result = append(result, testChr)
			// Remove from other list if present
			other = removeString(other, testChr)
		}
	}

	// Sort and add remaining chromosomes alphabetically
	sort.Strings(other)
	result = append(result, other...)

	return result
}

// removeString removes a string from a slice
func removeString(slice []string, s string) []string {
	for i, v := range slice {
		if v == s {
			return append(slice[:i], slice[i+1:]...)
		}
	}
	return slice
}

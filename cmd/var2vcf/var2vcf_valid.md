# var2vcf_valid — Go Implementation Status

Reference Perl: `AstraZeneca-NGS/VarDict` commit `009e017` (`var2vcf_valid.pl` v1.8.2)  
Go source: `cmd/var2vcf/` (main.go, header.go, sort.go)

---

## Status

| Component | Status | Notes |
|-----------|--------|-------|
| Column mapping (post-teststrandbias) | ✅ Done | Fixed critical bug where cols 20+ were all wrong |
| CLI flags & defaults | ✅ Done | All flags match Perl; defaults corrected |
| Depth filters (`d`, `v`) | ✅ Done | Uses `HiCnt` (col 30), not raw `VD` (col 8) |
| Quality filters (`f`, `p`, `q`, `Q`, `SN`, `NM`) | ✅ Done | All implemented with correct defaults |
| pSTD filter | ✅ Done | Default ON (1), matching Perl |
| Strand bias filter (`Bias`) | ✅ Done | All 5 Perl conditions implemented |
| MSI filter | ✅ Done | Non-monomer and monomer cases |
| LongMSI filter | ✅ Done | Both hiaf/msi threshold conditions |
| Cluster filter (`-c`) | ✅ Done | Deferred printing with retroactive marking of prev SNV |
| `-S` strict mode | ✅ Done | Skips non-PASS records |
| `-E` flag (omit END) | ✅ Done | Fixed inversion bug (was enabling, now disables END) |
| `-A` all-variants mode | ✅ Done | |
| Variant deduplication | ✅ Done | By `chr-start-end-ref-alt` key |
| AF-descending sort per position | ✅ Done | Matches Perl's `sort { $b->[14] <=> $a->[14] }` |
| SV handling (split/span reads) | ✅ Done | Reads from TAmp (col 37) |
| INFO fields (full set) | ✅ Done | Includes REFBIAS, VARBIAS, SHIFT3, MSI, MSILEN, HICNT, HICOV, LSEQ, RSEQ |
| SBF in INFO | ✅ Done | Fisher p-value float (was wrongly a strand counts string) |
| BIAS field transformation | ✅ Done | `;` → `:` in bias string |
| QUAL column computation | ✅ Done | `floor(log2(vd) * qual)` |
| FORMAT fields | ✅ Done | Reduced to Perl's 7 fields; RD/ALD as fwd,rev pairs |
| VCF header (INFO/FILTER/FORMAT defs) | ✅ Done | All definitions match Perl |
| Chromosome sort order | ✅ Done | Numeric → X/Y → MT/chrM → other |
| Amplicon mode (`-a`) | ⏳ TODO | GDAMP/TLAMP/NCAMP/AMPFLAG output stubbed; AMPBIAS filter not implemented |

---

## Known Limitations

### Amplicon Mode (TODO)
The `-a` flag and `AMPBIAS` filter are not implemented. Amplicon mode is detected
when the input has > 39 columns (col 39 = `$ampflag` is present). Currently:
- The `IsAmp` field on `VariantRecord` correctly detects amplicon input
- The GDAMP/TLAMP/NCAMP/AMPFLAG INFO fields are written when IsAmp=true
- The `AMPBIAS` filter logic (`gamp < tamp-ncamp || ampflag`) is **not** applied
- The `-a` flag does not exist yet

### Float Formatting
Go's `%g` is used for float output, which matches Perl's default scalar interpolation
for typical variant values. Edge cases with very high precision floats may produce
minor formatting differences that are numerically equivalent.

### Chromosome Sort Edge Case
Perl strips all non-digit characters to derive a numeric sort key (e.g. `GL000207.1`
→ key `2071`). Go uses `strconv.Atoi` after stripping the `chr` prefix, which fails
on non-standard contig names and falls back to alphabetical. Standard chromosomes
(1–22, X, Y, MT) sort identically.

---

## Files

| File | Purpose |
|------|---------|
| `main.go` | Core logic: parsing, filtering, VCF record construction |
| `header.go` | VCF header generation (INFO/FILTER/FORMAT definitions) |
| `sort.go` | Chromosome sort order |
| `Makefile` | Build targets |
| `var2vcf_valid_internals.md` | Detailed Perl vs Go gap analysis (reference doc) |
| `var2vcf_valid.md` | This file — implementation status |

---

## Pipeline Usage

```bash
vardict-java -G genome.fa -N sample -b regions.bed -f 0.01 input.bam \
  | teststrandbias \
  | var2vcf -N sample -f 0.01 > variants.vcf
```

## Building

```bash
cd cmd/var2vcf
make build
# or
go build -o var2vcf .
```

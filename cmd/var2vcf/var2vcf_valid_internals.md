# var2vcf_valid — Internals & Go vs Perl Gap Analysis

Reference Perl version: `AstraZeneca-NGS/VarDict` commit `009e017`  
Go version: `cmd/var2vcf/` in this repo

---

## Pipeline Context

`var2vcf_valid` is the **third stage** of the VarDict pipeline:

```
vardict-java  →  teststrandbias  →  var2vcf_valid
  (caller)         (strand bias)      (VCF output)
```

`teststrandbias` inserts two new columns (Fisher p-value and odds ratio) at positions 20–21 (0-indexed) into the raw VarDict output. `var2vcf_valid` reads this augmented stream.

---

## Input Column Format (post-teststrandbias)

The input is tab-separated with **36–42 columns** depending on VarDict mode. The canonical 40-column format is:

| Col (0-idx) | Variable | Description |
|-------------|----------|-------------|
| 0 | `sample` | Sample name |
| 1 | `gene` | Gene name |
| 2 | `chr` | Chromosome |
| 3 | `start` | Start position |
| 4 | `end` | End position |
| 5 | `ref` | Reference allele |
| 6 | `alt` | Alt allele |
| 7 | `dp` | Total depth |
| 8 | `vd` | Variant depth (raw read count) |
| 9 | `rfwd` | Ref forward reads |
| 10 | `rrev` | Ref reverse reads |
| 11 | `vfwd` | Variant forward reads |
| 12 | `vrev` | Variant reverse reads |
| 13 | `genotype` | Genotype string |
| 14 | `af` | Allele frequency |
| 15 | `bias` | Strand bias code (e.g. `2;1`) |
| 16 | `pmean` | Mean position in read |
| 17 | `pstd` | Position std dev in read |
| 18 | `qual` | Mean base quality |
| 19 | `qstd` | Base quality std dev |
| **20** | **`sbf`** | **Fisher p-value** ← inserted by teststrandbias |
| **21** | **`oddratio`** | **Odds ratio** ← inserted by teststrandbias |
| 22 | `mapq` | Mean mapping quality |
| 23 | `sn` | Signal to noise ratio |
| 24 | `hiaf` | High-quality allele frequency |
| 25 | `adjaf` | Adjusted AF (indel realignment) |
| 26 | `shift3` | Bases shifted 3' for deletions |
| 27 | `msi` | Microsatellite score |
| 28 | `msilen` | MSI unit length (bp) |
| 29 | `nm` | Mean mismatches in reads |
| 30 | `hicnt` | High-quality variant reads |
| 31 | `hicov` | High-quality total reads |
| 32 | `lseq` | 5' flanking sequence |
| 33 | `rseq` | 3' flanking sequence |
| 34 | `seg` | Segment (unused in var2vcf) |
| 35 | `type` | Variant type (SNV/Insertion/Deletion/Complex) |
| 36 | `gamp` | Supporting amplicons (or DUPRATE for non-amplicon) |
| 37 | `tamp` | Total amplicons (or SV split/span reads `SR-SP-cluster`) |
| 38 | `ncamp` | Non-working amplicons (or CRISPR field) |
| 39 | `ampflag` | Amplicon mismatch flag |

> **Note:** Presence of columns 36–39 determines amplicon mode. `defined($ampflag)` in Perl sets `$isamp = 1`.

---

## CLI Arguments

| Flag | Type | Perl default | Go default | Notes |
|------|------|-------------|------------|-------|
| `-N` | string | — | — | Sample name override |
| `-G` | string | — | — | Reference FASTA path (header only) |
| `-b` | string | — | — | BED file for contig headers |
| `-d` | int | **3** | 3 | Min total depth |
| `-v` | int | **2** | 2 | Min hi-qual variant depth |
| `-f` | float | **0.02** | 0.02 | Min allele frequency |
| `-p` | float | **8** | ~~5.0~~ | Min mean position in read |
| `-q` | float | **22.5** | 22.5 | Min mean base quality |
| `-Q` | float | **10** | ~~20.0~~ | Min mean mapping quality |
| `-o` | float | **1.5** | 1.5 | Min signal-to-noise ratio |
| `-F` | float | **0.2** | 0.2 | GT frequency threshold (hom/het) |
| `-P` | int (0/1) | **1** (ON) | ~~false (OFF)~~ | Filter pstd=0 variants |
| `-I` | int | **12** | missing | Max non-monomer MSI for AF<0.5 |
| `-m` | float | **5.25** | missing | Max mean mismatches |
| `-c` | int | **0** | missing | Cluster distance in bp |
| `-T` | int | **1** | 1 | Min split reads for SV |
| `-E` | bool | disables END tag | ~~enables END tag~~ | **Inverted in Go** |
| `-A` | bool | output all variants | output all variants | ✅ |
| `-S` | bool | strict (skip non-PASS) | missing | |
| `-a` | bool | amplicon mode | missing | |

---

## Processing Flow

### 1. Read all input into memory, group by `chr → position`

```
while(<STDIN>):
    skip lines containing "R_HOME"
    parse tab-separated fields
    group: hash{chr}{start} → list of variant arrays
    sort variants at each position by AF descending
```

### 2. Print VCF header

- `##fileformat=VCFv4.2`
- `##source=VarDict_v1.8.2`
- `##reference=` (if `-G` set)
- `##contig=` lines (if `-b` BED file set — reads chr + end from col 0,2)
- All `##INFO=` definitions
- All `##FILTER=` definitions (using runtime values of thresholds)
- All `##FORMAT=` definitions
- `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample>`

### 3. Process chromosomes in sorted order

Chromosome sort order:
1. Numeric chromosomes sorted numerically (chr1, chr2, …, chr22)
2. Sex chromosomes: X, Y (with or without `chr` prefix)
3. Mitochondrial: MT or chrM
4. All others alphabetically

> **Perl reorder quirk:** chromosomes with digits AND underscores (e.g. `chr1_random`) are treated as non-numeric and sorted alphabetically.

### 4. Process each position

For each position, variants are sorted by AF descending. Only the top variant is output unless `-A` is set. Deduplication by `chr-start-end-ref-alt` key prevents duplicate records.

### 5. Apply filters (build FILTER column)

Filters are **accumulated** (not short-circuit). A variant gets all applicable filter tags. If none apply, it gets `PASS`. With `-S`, non-PASS variants are skipped entirely.

The **Bias filter uses deferred printing**: the current variant is buffered as `$pinfo1/$pfilter/$pinfo2` and the *previous* variant is printed. This allows the Cluster filter to retroactively mark the previous SNV. The final variant is printed after the loop.

### 6. Build and print VCF record

```
CHROM  POS  .  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
```

---

## Filter Logic (Perl Reference)

### Depth filters
```perl
push @filters, "d$TotalDepth"  if $dp < $TotalDepth
    unless ($hicnt * $hiaf >= 0.5);    # exception for low-depth high-quality

push @filters, "v$VarDepth"    if $hicnt < $VarDepth    # uses hicnt (col 30), NOT vd (col 8)
    unless ($hicnt * $hiaf >= 0.5);
```

### Quality filters
```perl
push @filters, "f$Freq"   if $af < $Freq
push @filters, "p$Pmean"  if $pmean < $Pmean
push @filters, "pSTD"     if $opt_P && $pstd == 0 && !$isamp && $af < 0.35
push @filters, "q$qmean"  if $qual < $qmean
push @filters, "Q$Qmean"  if $mapq < $Qmean && $af < 0.8
push @filters, "SN$SN"    if $sn < $SN
push @filters, "NM$opt_m" if $nm > $opt_m
```

### MSI filters
```perl
# Non-monomer MSI
push @filters, "MSI$opt_I" if (
    ($msi > $opt_I && $msilen > 1 && $af < 0.2 && abs(len($ref)-len($alt)) == $msilen)
    || ($msi >= 13 && $msilen == 1 && $af <= 0.275 && abs(len($ref)-len($alt)) == $msilen)
)

# Long MSI (separate filter tag)
if abs(len($ref)-len($alt)) == $msilen:
    push @filters, "LongMSI" if $hiaf <= 0.275 && $msi >= 13
    push @filters, "LongMSI" if $hiaf <= 0.2   && $msi >= 8 && $msilen > 1
```

### Strand bias filter
```perl
push @filters, "Bias" if (
    $hiaf < 0.25
    && $bias eq "2;1"
    && $sbf < 0.01                       # Fisher p-value
    && ($oddratio > 5 || $oddratio == 0) # OR inverted to >= 1 before this check
    && ($end - $start) < 100
)
```

> `$oddratio` is pre-processed: `Inf → 0`, and if `0 < oddratio < 1` → `1/oddratio`.

### Cluster filter (SNV proximity)
```perl
# Applied only if no other filters and opt_c > 0 (default 0 = disabled)
if $type eq "SNV" && @filters == 0 && ($start - $pvs) < $opt_c:
    push @filters, "Cluster${opt_c}bp"
    # Also retroactively marks the previous SNV with same filter (via $pfilter)
```

### Amplicon bias filter
```perl
push @filters, "AMPBIAS" if $isamp && ($gamp < $tamp - $ncamp || $ampflag)
```

### SV skip (not a filter tag — silently skips the record)
```perl
next unless $splitreads >= $opt_T   # if alt contains '<' and split reads insufficient
```

---

## INFO Fields (Perl Reference)

```
SAMPLE=<sample_with_spaces_replaced_by_underscores>
TYPE=<SNV|Insertion|Deletion|Complex|REF>
DP=<dp>
END=<end>                           (omitted if -E flag set)
VD=<vd>
AF=<af>
BIAS=<bias with ; replaced by :>
REFBIAS=<rfwd>:<rrev>
VARBIAS=<vfwd>:<vrev>
PMEAN=<pmean>
PSTD=<pstd>
QUAL=<qual>
QSTD=<qstd>
SBF=<sbf>                           (Fisher p-value float, e.g. 0.03421)
ODDRATIO=<oddratio>
MQ=<mapq>
SN=<sn>
HIAF=<hiaf>
ADJAF=<adjaf>
SHIFT3=<shift3>
MSI=<msi>
MSILEN=<msilen>
NM=<nm>
HICNT=<hicnt>
HICOV=<hicov>
LSEQ=<lseq>
RSEQ=<rseq>
GDAMP=<gamp>  TLAMP=<tamp>  NCAMP=<ncamp>  AMPFLAG=<ampflag>  (amplicon mode only)
DUPRATE=<gamp>                      (non-amplicon mode, col 36)
SPLITREAD=<splitreads>  SPANPAIR=<spanpairs>  (if tamp is defined)
SVTYPE=<type>  SVLEN=<len>          (if alt contains '<')
```

> `END` is printed by **default**. The `-E` flag **disables** it.  
> `SBF` is the **p-value** (a float), not the strand counts string.

---

## FORMAT Fields (Perl Reference)

```
GT:DP:VD:AD:AF:RD:ALD
```

| Field | Value |
|-------|-------|
| `GT` | `0/0`, `0/1`, `1/0`, or `1/1` (see genotype logic below) |
| `DP` | `$dp` — total depth |
| `VD` | `$vd` — variant depth |
| `AD` | `$rd` if vd==0, else `$rd,$vd` (ref depth, then alt depth) |
| `RD` | `$rfwd,$rrev` — ref forward and reverse reads (comma-separated pair) |
| `ALD` | `$vfwd,$vrev` — alt forward and reverse reads (comma-separated pair) |

### Genotype logic
```perl
if ref == alt:
    alt = "."
    gt  = "0/0"
else:
    gt = (1-af) < GTFreq ? "1/1"
       : af >= 0.5       ? "1/0"
       : af >= Freq       ? "0/1"
       :                    "0/0"
```

### QUAL column
```perl
QUAL = (vd <= 1) ? 0 : int(log2(vd) * qual)
```

---

## Key Bugs in the Current Go Implementation

| # | Severity | Issue |
|---|----------|-------|
| 1 | **Critical** | Column mapping wrong from col 20 onwards — entire MapQ, SN, HiAF, PValue etc. read from wrong positions |
| 2 | **Critical** | `-E` flag is **inverted**: Perl `-E` disables END tag; Go `-E` enables it |
| 3 | **Critical** | `SBF` in INFO outputs strand counts string instead of Fisher p-value float |
| 4 | **High** | Bias filter logic is completely different (Go: `pvalue<0.01 && OR>=2`; Perl has 5 conditions) |
| 5 | **High** | Variant depth filter uses `AltDepth` (col 8 = raw vd); should use `hicnt` (col 30 = hi-qual reads) |
| 6 | **High** | `pSTD` filter default is inverted (Go: OFF; Perl: ON) |
| 7 | **High** | Default `-p` is 5.0; should be **8** |
| 8 | **High** | Default `-Q` is 20.0; should be **10** |
| 9 | **High** | Missing `NM` filter (`nm > 5.25`) |
| 10 | **High** | Missing `MSI` and `LongMSI` filters |
| 11 | **Medium** | Missing INFO fields: `REFBIAS`, `VARBIAS`, `SHIFT3`, `MSI`, `MSILEN`, `HICNT`, `HICOV`, `LSEQ`, `RSEQ` |
| 12 | **Medium** | `BIAS` field in INFO needs `;` → `:` transformation |
| 13 | **Medium** | FORMAT has wrong fields — Perl uses 7; Go has 20 with wrong `RD`/`ALD` format |
| 14 | **Medium** | `QUAL` column is always `.`; should be computed as `log2(vd) * qual` |
| 15 | **Medium** | Missing `-S` strict mode (skip non-PASS variants) |
| 16 | **Medium** | No variant deduplication by `chr-start-end-ref-alt` |
| 17 | **Low** | Missing `Cluster` filter (default opt_c=0 means disabled, so low impact) |
| 18 | **Low** | Missing `-a` amplicon mode and `AMPBIAS` filter |
| 19 | **Low** | Missing `-I`, `-m`, `-c` CLI flags |
| 20 | **Low** | SV split/span reads read from wrong column |

# var2vcf - Go Implementation

A Go implementation of VarDict's `var2vcf_valid.pl` (v1.8.2), converting VarDict variant output to VCF format. Drop-in replacement for the Perl script in non-amplicon workflows.

## Quick Start

```bash
# Build
make build

# Use in pipeline
vardict-java -G genome.fa -N sample -b regions.bed -f 0.01 input.bam \
  | teststrandbias \
  | ./var2vcf -N sample -f 0.01 > variants.vcf
```

## Command-Line Options

```
-N string   Sample name override
-G string   Path to reference FASTA (written to ##reference header)
-b string   Path to BED file (written as ##contig headers)

-d int      Minimum total depth (default: 3)
-v int      Minimum high-quality variant depth (default: 2)
-f float    Minimum allele frequency (default: 0.02)
-p float    Minimum mean position in read (default: 8)
-q float    Minimum mean base quality (default: 22.5)
-Q float    Minimum mean mapping quality (default: 10)
-o float    Minimum signal-to-noise ratio (default: 1.5)
-F float    Genotype frequency threshold for homozygous call (default: 0.2)
-I int      Maximum non-monomer MSI for AF<0.5 variants (default: 12)
-m float    Maximum mean mismatches in reads (default: 5.25)
-c int      Filter SNVs within this many bp of each other, 0=disabled (default: 0)
-P int      Filter variants with pstd=0: 1=yes, 0=no (default: 1)
-T int      Minimum split reads for structural variants (default: 1)

-E          Do not print END tag in INFO field
-A          Output all variants at the same position (default: highest AF only)
-S          Strict mode: only output variants that pass all filters
```

## Filters Applied

| Filter tag | Condition |
|-----------|-----------|
| `d<N>` | Total depth < `-d` |
| `v<N>` | High-quality variant depth < `-v` |
| `f<N>` | Allele frequency < `-f` |
| `p<N>` | Mean position in read < `-p` |
| `pSTD` | Position std dev = 0 (when `-P 1`) |
| `q<N>` | Mean base quality < `-q` |
| `Q<N>` | Mean mapping quality < `-Q` and AF < 0.8 |
| `SN<N>` | Signal-to-noise < `-o` |
| `NM<N>` | Mean mismatches > `-m` |
| `MSI<N>` | Variant in microsatellite region |
| `LongMSI` | Variant flanked by long MSI run |
| `Bias` | Strand bias (hiaf<0.25, bias=2;1, p<0.01, OR>5 or OR=0, len<100) |
| `Cluster<N>bp` | Two SNVs within `-c` bp of each other |

## Limitations

- **Amplicon mode** (`-a` flag) is not yet implemented. AMPBIAS filter and amplicon-specific INFO fields (GDAMP/TLAMP/NCAMP/AMPFLAG) are pending.
- Designed for the standard single-sample (`var2vcf_valid.pl`) workflow, not the paired tumor-normal workflow (`var2vcf_paired.pl`).

## License

MIT License (same as VarDict)

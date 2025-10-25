# var2vcf - Go Implementation

A fast, memory-efficient Go implementation of VarDict's `var2vcf_valid.pl` script for converting variant calls to VCF format.

## Quick Start

```bash
# Build
make build

# Run tests
make test

# Use in pipeline
cat input.tsv | ./var2vcf -N sample_name -f 0.01 > output.vcf
```

## Features

✅ **Zero runtime dependencies** - Static binary, no Perl required  
✅ **Memory efficient** - Processes chromosomes one at a time  
✅ **Fast execution** - Go's performance advantages  
✅ **Compatible** - Drop-in replacement for `var2vcf_valid.pl`  
✅ **Well-tested** - Comprehensive unit tests  

## Usage

```bash
# Basic usage
./var2vcf -N sample_name -f 0.01 < input.tsv > output.vcf

# With full options
./var2vcf \
  -N MySample \
  -G /path/to/genome.fa \
  -b /path/to/regions.bed \
  -f 0.01 \
  -d 8 \
  -v 4 \
  -E \
  < input.tsv > output.vcf
```

## Command-Line Options

```
-N string   Sample name (required)
-G string   Reference genome path
-b string   BED file path for contig information
-d int      Minimum total depth (default: 8)
-v int      Minimum variant depth (default: 4)
-f float    Minimum allele frequency (default: 0.01)
-p float    Minimum mean position in read (default: 5.0)
-q float    Minimum base quality (default: 22.5)
-Q float    Minimum mapping quality (default: 20.0)
-S float    Minimum signal-to-noise ratio (default: 1.5)
-F float    Genotype frequency threshold (default: 0.2)
-E          Print END tag in INFO field
-A          Print all variants at same position
-P          Filter variants with position std dev = 0
-T int      Minimum split reads for SVs (default: 1)
```

## Complete VarDict Pipeline

```bash
# Example pipeline
vardict-java -G genome.fa -N sample -b input.bam -f 0.01 regions.bed \
  | ./teststrandbias \
  | ./var2vcf -N sample -f 0.01 > variants.vcf
```

## Documentation

See the full documentation at [https://github.com/yourusername/var2vcf](https://github.com/yourusername/var2vcf)

## License

MIT License (same as VarDict)

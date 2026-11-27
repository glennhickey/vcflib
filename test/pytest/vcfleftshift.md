% VCFLEFTSHIFT(1) vcfleftshift | Trim and left-align variants without moving them across a haplotype's other variants
% Glenn Hickey

# NAME

**vcfleftshift** trims and left-aligns variants against the reference, like `bcftools norm -f`, but never moves a variant onto or across a base that another variant on one of its haplotypes describes.

# SYNOPSIS

**vcfleftshift** -f *ref.fa* < *in.vcf* > *out.vcf*

# DESCRIPTION

`bcftools norm -f` places each record by the reference alone. In a phased VCF that can slide an insertion onto the base of a SNP the same haplotype carries, or a deletion through a tandem copy and across that haplotype's own SNPs. That haplotype then has two records describing one base, or a SNP inside its own deletion.

**vcfleftshift** keeps, for each haplotype, the end of the last variant it carries, and moves an indel left only while it stays to the right of that end for every haplotype that carries it. The bases it moves across are reference on all of those haplotypes, so the move changes none of their sequences. Records at the same position that share a haplotype are trimmed but not moved. Multiallelic and symbolic records are passed through: split them with `bcftools norm -m -any` first. The input must be sorted by position, and moved records can come out of order, so sort the output afterwards.

## Options

```python
>>> from rtest import run_stdout, head, cat, sh
>>> head("vcfleftshift -h", 4)
usage: vcfleftshift -f ref.fa [options] < in.vcf > out.vcf
>
Trim and left-align variants against the reference, like bcftools norm -f, but never
move a variant onto or across a base that another variant on one of its haplotypes

```

## An insertion next to a SNP

The reference has `GGAGAGAGA` at 11-19. S1's first haplotype has a G>A SNP at 11 and inserts `GA` into the repeat. S2's first haplotype inserts `GAGA` there and has no SNP. Both insertions are written at the right end of the repeat, as vcfwave places them. `bcftools norm -f` would move both to 11, putting S1's insertion on the base its SNP changes. vcfleftshift stops S1's insertion at 12, after the SNP, and moves S2's all the way:

```python
>>> sh("vcfleftshift -f ../samples/vcfleftshift.fa < ../samples/vcfleftshift.vcf | grep -v '^#' | cut -f1-5,10,11 | tr '\\t' ' '")
chr1 11 snp G A 1|0 0|0
chr1 12 ins_GA G GAG 1|0 0|0
chr1 11 ins_GAGA G GGAGA 0|0 1|0

```

Without genotypes nothing limits a move, so both insertions go to 11, as with `bcftools norm -f`:

```python
>>> sh("cut -f1-8 ../samples/vcfleftshift.vcf | vcfleftshift -f ../samples/vcfleftshift.fa | grep -v '^#' | cut -f1-5 | tr '\\t' ' '")
chr1 11 snp G A
chr1 11 ins_GA G GGA
chr1 11 ins_GAGA G GGAGA

```

# LICENSE

Copyright 2026 Glenn Hickey and vcflib contributors. MIT licensed.

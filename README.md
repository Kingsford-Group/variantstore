# VariantStore

VariantStore: A Large-Scale Genomic Variant Search Index

Overview
--------

VariantStore is a system for efficiently indexing and querying genomic
information (genomic variants and phasing information) from thousands of samples
containing millions of variants. VariantStore builds a variant-oriented index by
mapping each variant to the list of samples that contain the variant.  Variants
are indexed based on the positions where they occur. It supports querying
variants occurring between two positions across a chromosome based on the
reference or a sample coordinate. VariantStore can scale to tens of millions of
variants and thousands of samples and can efficiently scale out-of-RAM to
storage devices to enable memory-efficient construction and query.

API
--------
* `variantstore construct`: construct a variation graph and position index.
* `variantstore query`: query variation graph for variants using the position index.

Build
-------
Library dependencies (given version or higher):
- [protobuf 3.7.1](https://github.com/protocolbuffers/protobuf)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)

```bash
$ make proto
$ make variantstore
```

Construct variation graph

```bash
$ ./variantstore construct -r data/x.fa -v data/x.vcf.gz -p ser/
```

The usage for this command are as follows:

```bash
SYNOPSIS
        variantstore construct -r <reference-file> -v <vcf-file> -p <output-prefix>

OPTIONS
        <reference-file>
                    reference file

        <vcf-file>  variant call format file

        <output-prefix>
                    output directory
```

Variant queries

```bash
$ ./variantstore query -p ser/ -t 6 -r 10:105 -m 0 -s 1 -o variant.txt -v
```

The usage for this command are as follows:

```bash
SYNOPSIS
        variantstore query -p <output-prefix> -t <query-type> -r <region> -m <mode> [-o <outfile>] [-s <sample-name>] [-a <alt-seq>] [-b <ref-seq>] [-v]

OPTIONS
        <output-prefix>
                    output directory

        <query-type>
                    Types of query.
                   1. Get sample's sequence in reference coordinates.
                   2. Get sample's sequence in sample coordinates.
                   3. Return closest variant in reference coordinates.
                   4. Get sample's variants in reference coordinates.
                   5. Get sample's variants in sample coordinates.
                   6. Get variants in reference coordinates.
                   7. Return samples with a given mutation.

        <region>    region in format <start>:<end>, regions separated by ','
        <mode>      READ_INDEX_ONLY: 0, READ_COMPLETE_GRAPH:1
        <outfile>   output_file

        <sample-name>
                    sample name

        <alt-seq>   alternative sequences, separated by ','
        <ref-seq>   reference sequences, separated by ','

        -v, --verbose
                    print vcf
```

Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey2@cs.cmu.edu>
- Yinjie Gao <yinjieg@andrew.cmu.edu>
- Carl Kingsford <carlk@cs.cmu.edu>



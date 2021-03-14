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
- g++ 7.5.0 (VariantStore needs C++17. C++17 features are available since GCC 5. But not tested with GCC version < 7.5.0.)
- [protobuf 3.7.1](https://github.com/protocolbuffers/protobuf)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- zlib1g-dev 1:1.2.11.dfsg-0ubuntu2
- libbz2-dev 1.0.6-8.1ubuntu0.2
- liblzma-dev 5.2.2-1.3
```bash
$ make proto
$ make variantstore
```

Construct variation graph

```bash
$ ./variantstore construct -r data/x.fa -v data/x.vcf.gz -p ser/
```

Expected outputv
```bash
[2020-05-19 10:10:35.172] [info] Creating variant graph
[2020-05-19 10:10:35.208] [info] Building sample vector based variant graph.
[2020-05-19 10:10:35.208] [info] Adding mutations from: data/x.vcf.gz #Samples: 1
[2020-05-19 10:1u0:35.211] [info] Num mutations: 75 num mutations-sample: 75
[2020-05-19 10:10:35.211] [info] Num vars: 75
[2020-05-19 10:10:35.212] [info] Fixing sample indexes in the graph
[2020-05-19 10:10:35.212] [info] Graph stats:
[2020-05-19 10:10:35.212] [info] Chromosome: x #Vertices: 212 #Edges: 287 Seq length: 1074
[2020-05-19 10:10:35.212] [info] Serializing variant graph to disk
[2020-05-19 10:10:35.245] [info] Number of sample vector classes: 2
[2020-05-19 10:10:35.245] [info] Creating Index
[2020-05-19 10:10:35.246] [info] Serializing index to disk
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

Expected output
```bash
[2020-05-19 10:19:08.284] [info] Loading Index ...
[2020-05-19 10:19:08.284] [info] Loading variant graph ...
[2020-05-19 10:19:08.284] [info] Read index only ..
[2020-05-19 10:19:08.329] [info] Graph stats:
[2020-05-19 10:19:08.329] [info] Chromosome: x #Vertices: 212 #Edges: 0 Seq length: 1074
[2020-05-19 10:19:08.329] [info] 6. Get variants in ref coordinate. 0
Number of variants get_var_in_ref: 8
Query1: (query_var_in_ref) Total Time Elapsed: 0.000426seconds
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

Visualize variation graph

Output variation graph in the dot format

```bash
$ ./variantstore draw -p ser/ -r 10 -h 3 
```

Visualize dot graph

```bash
$ make graph
```

Expected output
```bash
[2021-03-14 15:23:18.632] [info] Loading Index ...
[2021-03-14 15:23:18.632] [info] Loading variant graph ...
[2021-03-14 15:23:18.632] [info] Read complete graph ..
[2021-03-14 15:23:18.682] [info] Graph stats:
[2021-03-14 15:23:18.682] [info] Chromosome: x #Vertices: 212 #Edges: 0 Seq length: 1074
[2021-03-14 15:23:18.682] [info] Looking up vertex corresponding to the queried region
```

The usage for this command are as follows:

```bash
SYNOPSIS
        variantstore draw -p <output-prefix> -r <region> -h <hops> [-s <sample-name>] [-o <outfile>]

OPTIONS
        <output-prefix>
                    output directory

        <region>    region in format <start>:<end>, regions separated by ','
        <hops>      radius of the subgraph

        <sample-name>
                    sample name (default: REF)

        <outfile>   output_file

```



Contributing
------------
Contributions via GitHub pull requests are welcome.


Authors
-------
- Prashant Pandey <ppandey2@cs.cmu.edu>
- Yinjie Gao <yinjieg@andrew.cmu.edu>
- Carl Kingsford <carlk@cs.cmu.edu>



# What is it?

The CracTools-SimCT software (or simply SimCT) is a complete and modular workflow to simulate RNA-Seq data from a reference genome and known transcript annotations.

SimCT works in three steps :

1. The first step introduces a set of variants in the reference genome. SimCT then generates a haploid mutated genome in FASTA format and a set of GTF annotations whose coordinates are converted for this modified reference.
2. We then process the modified reference files using FluxSimulator, that offers good performance in generating a random expression profile and simulating reads by reproducing a complete ”in silico” RNA-Seq protocol.
3. In the final step, the reads and corresponding alignments produced by FluxSimulator are sent to post-processing where alignment coordinates are converted to the coordinates of the original reference genome. The errors are then extracted from read sequence (encoded in lowercase by FluxSimulator) and a new FASTQ file is produced with alignments and errors encoded in the read name.

# Table of contents
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [What is it?](#what-is-it)
- [Table of contents](#table-of-contents)
- [Installation](#installation)
	- [Requirements](#requirements)
	- [Install from tarball](#install-from-tarball)
	- [Install from sources](#install-from-sources)
- [Quick usage examples](#quick-usage-examples)
- [Documentation](#documentation)
	- [Input files](#input-files)
		- [Genome directory](#genome-directory)
		- [Annotations](#annotations)
	- [Output files](#output-files)
	- [Options](#options)
	- [Read name encoding](#read-name-encoding)
- [Authors](#authors)

<!-- /TOC -->

# Installation

## Requirements

* Perl5 distribution
* [cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7043/lib/App/cpanminus.pm) (CPAN minus)
* [FluxSimulator]('http://artifactory.sammeth.net/artifactory/barna/barna/barna.simulator/1.2.1/flux-simulator-1.2.1.tgz') Make sure that the binary "flux-simulator" is available in the $PATH to be executed by SimCT.
* [CracTools-core](https://metacpan.org/release/CracTools) perl package. It will be automaically installed by cpanm along with all other CPAN dependancies.

## Install from tarball

This is the simpliest way to install simCT.

1. Go to the [release page](https://github.com/jaudoux/simct/releases) of the github projet.
1. Download the latest tarball (tar.gz) release : `wget https://github.com/jaudoux/simct/releases/download/$VERSION/CracTools-SimCT-$VERSION.tar.gz`
2. Install the package with cpanm : `cpanm [-l local_dir] CracTools-SimCT.tar.gz`

If you do not have admin rights, you can use the option `-l` to specify cpanm a local directory to install simCT.

## Install from sources

To install SimCT from the sources, you will need [Dist::Zilla]('http://dzil.org/') software, which is a Perl package manager.

1. Clone SimCT repository : `git clone https://github.com/jaudoux/simct.git`
2. Build and install : `dzil install --install-command 'cpanm [-l LOCAL_INSTALLATION_DIRECTORY].'` (ommit the `-l` option if you want to install simCT to the system).

# Quick usage examples


**Basic 1M Paired-end reads simulation**:<br>
```simCT -g genome_dir -a annotations.gtf --nb-reads 1000000```

**40M reads with 100 fusion genes**:<br>
`simCT -g genome_dir -a annotations.gtf --nb-reads 40000000 --nb-fusions 100`

**40M reads with 95% of mutations taken from a VCF files**:<br>
`simCT -g genome_dir -a annotations.gtf --nb-reads 40000000 `

# Documentation

In order to simulate RNA-Seq datasets with SimCT, all you need is a reference genome as a set a FASTA files for each chromosome and a GTF file containing annotations of transcripts. Then run simCT as follow by indicationg the reference genome directory (where the FASTA are located), the GTF file, and the output directory :

    simCT -g reference_genome/ -a annotations.gtf [-o my_simulation]

[Many options](#options) are available, run `simCT --help` for more documentation.

## Input files
The mandatory input files for simCT are a directory containaing the reference genome as a set of FASTA files for each chromosome and a GTF file containing the annotations.

### Genome directory

SimCT expect the genome to be a directory that contains a FASTA file for each chromosome of the genome. For example :


    ├── genome_dir
    │   ├── chr1.fa
    │   ├── chr2.fa
    │   ├── chr3.fa

If you have a genome in multi-FASTA format, you can use this [splitFasta.pl](scripts/splitFasta.pl) perl script] to split it in multiple FASTA files.


    # Create a directory for your genome_dir
    mkdir genome_dir && cd genome_dir

    # Invoke the splitFasta script
    ./splitFasta.pl genome.fa


### Annotations

SimCT should work with all GTF files (Ensembl, UCSC, GENCODE, ...), since we use only generic fields (*chr*, *start*, *end*, *strand*, *gene_id*, *transcript_id*) and only exon features.

## Output files

The output directory of a simCT simulation should be structured :

```
├── gene-counts.tsv.gz
├── transcript-counts.tsv.gz
├── mutations.vcf.gz
├── reads_1.fastq.gz
├── reads_2.fastq.gz
├── splices.bed.gz
├── chimeras.tsv.gz
├── info.txt
├── FluxSimulator/
├── simulated_genome/
│   ├── chr{N}.fa
│   ├── annotations.gtf
│   ├── fusion{N}.fa
```

The following table describes the content of each of these files/directories.

File name | Description
------------|-------------
gene-counts.tsv.gz | A tabulated file with two colomns : `feature` and `truth`. Feature old the `gene_id` from the GTF add truth the number of read mapped to this gene.
transcript-counts.tsv.gz | A tabulated file with two colomns : `feature` and `truth`. Feature old the `transcript_id` from the GTF add truth the number of read mapped to this transcript.
mutations.vcf.gz | A [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) containing the inserted mutations that are supported by 1 read at least. The `ID` column contains the IDs of the reads supporting the mutation  (separated with a colon).
splices.bed.gz | A 6 column bed file containing the splice junctions supported by at least 1 read. The `name` column contains the IDs of the reads supporting the splice junction (separated with a colon).
chimeras.gz | A tabulated file (with no header line) containing the chimeric junctions. The file contains 7 columns. The first 6 column describes the breakpoint : `chr1`, `pos1`, `strand1`, `chr2`, `pos2`, `strand2`. The last column contains the IDs of the reads supporting the chimeric junction (separated with a colon).
reads_{1,2}.fastq.gz | The FASTQ file(s) holding the generated reads. Read names holds the truth for read alignments to the reference genome. Read the '[read name encoding](#read-name-encoding)'' section for more informations.
info.txt | A text file describing the number of element for each feature (read, SNV, Splice, Fusion, SNV, Indel) generated in the simulations.
**FluxSimulator/** | The directory contains [output files](http://sammeth.net/confluence/display/SIM/Appendix+A+-+File+Format+Specifications) generated by FluxSimulator.
**simulated_genome/** | The directory contains the fasta sequences of the simulated genome as well as the annotations that have been liftover this genome. If Fusions have been generated, each fusion will be printed in a separated fasta file.

## Options

Option name | Default value | Description
------------|---------------|-------------
-o,--output-dir | `simCT_simulation` | Location of the outout directory
-g,--genome-dir | **mandatory** | xReference genome directory (with chromosomes splited in individual FASTA files)
-a,--annotations | **mandatory**  | Annotations file in GTF format
-s,--substitution-rate | 0.001 | Rate a wich substitutions are randomly inserted in the reference genome
-i,--insertion-rate     | 0.0005 | Rate a wich insertions are randomly inserted in the reference genome
-d,--deletion-rate      | 0.0005 | Rate a wich deletions are randomly inserted in the reference genome
-f,--nb-fusions         | 0 | Number of fusions introduced in the simulated genome
--vcf-file              | NA | A VCF file that contains mutations
--vcf-ratio             | NA | Ratio of mutations that taken from the VCF source instead of random(default: 0.8)
--flux-par              | NA | Flux parameter file (prioritary over the following option)
--single-end            | NA| Single-end protocol
--nb-molecules          | 5000000 | Number of moleclules in the sample
--nb-reads              | 1000000| Number of reads sequenced
--reads-length          | 100 | Length of sequenced reads
--fragment-length       | 300 | Mean fragment length
--fragment-sd           | 75 | Standard deviation of fragment length
--disable-error-encoding | NA | Remove error encoding from read names
--uniq-ids               | NA | Read names are identical for both pairs

## Read name encoding
SimCT encodes the read alignement and error position inside the read name. This way, it is very easy to verify an alignment from a SAM entry.

The read format is defined as : `read_id:(chr,(-)pos,cigar(;)?)+,base64(err_pos?)+`

The read name is composed of three components:

1. the read id from 0 to `nb_reads - 1`
2. a list of SAM-like read alignments (chr, starting position and cigar chain)
3. the positions of sequencing errors (if any) encoded in base64.

For paired-end reads, a single read name is generated for both reads by using a single `read\_id`, concatenating the alignments and merging the positions of sequencing errors.


Example of a read name: `0:1,12345,20M2I3X30M;2:-6789:20M:CABgAAAKAAEAABgCAAAABABAQAAgAAACDD`

The encoded error list (`CABgAAAKAAEAABgCAAAABABAQAAgAAACDD`) is decoded as as follow : 1,12,23,43,45,62,78,89,91,120,132,148,167,187,192,193,198,199

# Authors
Jérôme Audoux - jaudoux@cpan.org
Nicolas Philippe - nphilippe@cpan.org
Mikaël Salson - mikael.salson@univ-lille1.fr

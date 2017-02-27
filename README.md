What is it?
==========
The CracTools-SimCT software (or simply SimCT) is a complete and modular workflow to simulate RNA-Seq data from a reference genome and known transcript annotations.

SimCT works in three steps :

1. The first step introduces a set of variants in the reference genome. SimCT then generates a haploid mutated genome in FASTA format and a set of GTF annotations whose coordinates are converted for this modified reference.
2. We then process the modified reference files using FluxSimulator, that offers good performance in generating a random expression profile and simulating reads by reproducing a complete ”in silico” RNA-Seq protocol.
3. In the final step, the reads and corresponding alignments produced by FluxSimulator are sent to post-processing where alignment coordinates are converted to the coordinates of the original reference genome. The errors are then extracted from read sequence (encoded in lowercase by FluxSimulator) and a new FASTQ file is produced with alignments and errors encoded in the read name.


Installation
============
Installing simCT is very simple, download the tarball and install SimCT
with the following command :

    cpanm CracTools-SimCT.tar.gz

If you do not have admin rights, you can use the option `-l` to specify cpanm a
local directory to install benchCT.

SimCT has a strong dependance to CracTools-core distribution wich is freely
available at CPAN. If you use `cpanm` they will be automatically installed.

You also need to install [FluxSimulator]('http://artifactory.sammeth.net/artifactory/barna/barna/barna.simulator/1.2.1/flux-simulator-1.2.1.tgz'), and make sure that the binary "flux-simulator" is available in the $PATH to be executed by SimCT.

Documentation
=============

In order to simulate RNA-Seq datasets with SimCT, all you need is a reference genome as a set a FASTA files for each chromosome and a GTF file containing annotations of transcripts. Then run simCT as follow by indicationg the reference genome directory (where the FASTA are located), the GTF file, and the output directory :

    simCT -g reference_genome/ -a annotations.gtf [-o my_simulation]

Many options are available, run `simCT --help` for more documentation.

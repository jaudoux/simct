use strict;
use warnings;

use Test::More tests => 1;
use CracTools::Utils;
use CracTools::SimCT::FluxWrapper;
use Inline::Files 0.68;
use File::Temp;

{
  # Load the fasta file
  my $genome_dir = File::Temp->newdir();
  my $chr1_fasta_file = CracTools::Utils::getWritingFileHandle("$genome_dir/1.fa");
  while(<CHR1_FASTA>) {print $chr1_fasta_file $_;}
  close $chr1_fasta_file;

  # Load the GTF annotations
  my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<ANNOTATIONS>) {print $gtf_file $_;}
  close $gtf_file;

  #my $fw = CracTools::SimCT::FluxWrapper->new();

  #my $output_dir = File::Temp->newdir();
  #$fw->generateSimulation(
  #  genome_dir      => $genome_dir,
  #  annotation_file => $gtf_file,
  #  output_dir      => $output_dir,
  #  flux_parameters => {
  #    POLYA_SHAPE   => 'NaN',
  #    TSS_MEAN      => 'NaN',
  #    READ_NUMBER   => 1000,
  #    #READ_LENGTH   => 75,
  #    PAIRED_END    => 'YES',
  #    NB_MOLECULES  => 100,
  #  },
  #);
  ok(1);
}

# Exon 1: CCCGTCGCATG
# Exon 2; AAGTCGAATGA

__CHR1_FASTA__
>1
ATGGTAGTACCCGTCGCATGTCGAAAGTGTCGAAAGTCCAAGTCGAATGATGGTAGTA
__ANNOTATIONS__
1	StringTie	transcript	10	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	10	20	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	40	50	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "2"; cov "161.272552";

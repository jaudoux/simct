use strict;
use warnings;

use Test::More tests => 9;
use CracTools::SimCT::GenomeSimulator;
use CracTools::Utils;
use Inline::Files 0.68;
use File::Temp;

# Load the fasta file
my $chr1_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
while(<CHR1_FASTA>) {print $chr1_fasta_file $_;}
close $chr1_fasta_file;
my $chr2_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
while(<CHR2_FASTA>) {print $chr2_fasta_file $_;}
close $chr2_fasta_file;

# Load the GTF annotations
my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
while(<ANNOTATIONS>) {print $gtf_file $_;}
close $gtf_file;

# Create the GenomeSimulator
my $gs = CracTools::SimCT::GenomeSimulator->new(
  reference_sequence_files => { 1 => $chr1_fasta_file->filename, 2 => $chr2_fasta_file->filename},
  annotation_file => $gtf_file,
);

$gs->addFusion('geneA','11,27','geneB','16,22');

my $genome_dir = File::Temp->newdir();
my $gtf_output = new File::Temp( SUFFIX => '.gtf' );

my $frozen_gs = $gs->generateGenome(
  fasta_dir   => $genome_dir,
  gtf_file    => $gtf_output,
  #gtf_file    => "test.gtf",
);

my $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/chrFusions.fa");
my $entry     = $fasta_it->();


# CTAGCTAGTTAGCTCGATCGGC => GCCGATCGAGCTAACTAGCTAG
my $fusion_1 = "(TGGTAGTACCCGTCGCATGTCGAAAGT)(GCCGATCGAGCTAACTAGCTAG)";
ok($entry->{seq} =~ /^$fusion_1$/,"generateGenome - FASTA control (2)");

my $gtf_it        = CracTools::Utils::gffFileIterator($gtf_output,'gtf');

$gtf_it->() for 1..5;

# Check the first fusion
{
  my $fusion_exon_1   = $gtf_it->(); # exon1, geneA
  my $fusion_exon_2   = $gtf_it->(); # exon2, geneA
  my $fusion_exon_3   = $gtf_it->(); # Fused exon (exon3, geneA) <=> (exon2, geneB)
  my $fusion_exon_4   = $gtf_it->(); # exon1, geneB

  is($fusion_exon_1->{start},1);
  is($fusion_exon_1->{end},6);
  is($fusion_exon_2->{start},8);
  is($fusion_exon_2->{end},9);
  is($fusion_exon_3->{start},11);
  is($fusion_exon_3->{end},34);
  is($fusion_exon_4->{start},39);
  is($fusion_exon_4->{end},49);
}

__CHR1_FASTA__
>chr1
ATGGTAGTACCCGTCGCATGTCGAAAGTCT
__CHR2_FASTA__
>chr2
GCTAGCTAGTTAGCTCGATCGGCAT
__ANNOTATIONS__
1	StringTie	exon	2	7	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	9	10	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	12	28	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	2	12	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	17	23	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";

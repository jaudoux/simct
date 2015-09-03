use strict;
use warnings;

use Test::More tests => 5;
use CracTools::SimCT::GenomeSimulator;
use CracTools::Utils;
use Inline::Files 0.68;
use File::Temp;

{
  # Load the fasta file
  my $fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
  while(<CHR1_FASTA>) {print $fasta_file $_;}
  close $fasta_file;

  # Load the GTF annotations
  my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<ANNOTATIONS>) {print $gtf_file $_;}
  close $gtf_file;

  # Create the GenomeSimulator
  my $gs = CracTools::SimCT::GenomeSimulator->new(
    reference_sequence_files => { 1 => $fasta_file },
    annotation_file => $gtf_file,
  );

  is($gs->getReferenceLength(1),30);

  # Add mutations
  $gs->addSubstitution("1",10,"G");
  $gs->addDeletion("1",15,10);
  $gs->addInsertion("1",2,"AGG");

  my $genome_dir = File::Temp->newdir();
  my $gtf_output = new File::Temp( SUFFIX => '.gtf' );

  $gs->generateGenome(
    fasta_dir   => $genome_dir,
    gtf_file    => $gtf_output,
    #gtf_file    => "test.gtf",
  );

  # Verify if fasta is good
  my $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/chr1.fa");
  my $entry     = $fasta_it->();
  is($entry->{seq},"ATAGGGGTAGTACGCGTCAGTCT","generateGenome - FASTA control (1)");

  # Verify if annotations are good
  my $gtf_it        = CracTools::Utils::gffFileIterator($gtf_output,'gtf');
  my $first_exon    = $gtf_it->();
  my $second_exon   = $gtf_it->();
  is($first_exon->{start},13,"generateGenome - GTF control (1)");
  is($first_exon->{end},18,"generateGenome - GTF control (2)");
  # Exon 2 is supressed by the deletion
  is($second_exon, undef,"generateGenome - GTF control (3)"); 
}

__CHR1_FASTA__
>chr1
ATGGTAGTACCCGTCGCATGTCGAA
AGTCT
>chr2
GCTAGCTAGTTAGCTCGATCGGCAT
__ANNOTATIONS__
1	StringTie	transcript	10	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	10	15	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	20	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "161.272552";
use strict;
use warnings;

use Test::More tests => 6;
use CracTools::SimCT::Genome;
use CracTools::Utils;
use Inline::Files 0.68;
use File::Temp;

{
  # Load the fasta file
  my $chr1_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
  while(<CHR1_FASTA>) {print $chr1_fasta_file $_;}
  close $chr1_fasta_file;
  my $chr2_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
  while(<CHR2_FASTA>) {print $chr2_fasta_file $_;}
  close $chr2_fasta_file;
  
  my $genome = CracTools::SimCT::Genome->new(
    reference_sequence_files => { "1" => $chr1_fasta_file->filename, "2" => $chr2_fasta_file->filename} ,
  );

  is($genome->getReferenceLength(1),30);
  is($genome->getReferenceLength(2),20);
  is($genome->getReferenceFile(1),$chr1_fasta_file->filename);
  is($genome->getReferenceFile(2),$chr2_fasta_file->filename);
  is($genome->getReferenceSeq(1),"ATGGTAGTACCCGTCGCATGTCGAAAGTCT");

  is(scalar $genome->references,2);
}

__CHR1_FASTA__
>chr1
ATGGTAGTACCCGTCGCATGTCGAA
AGTCT
__CHR2_FASTA__
>chr2
GCTAGCTAGTTAGCTCGATC

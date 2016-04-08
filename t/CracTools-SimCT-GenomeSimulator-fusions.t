use strict;
use warnings;

use Test::More tests => 9;
use CracTools::SimCT::Const;
use CracTools::SimCT::Genome;
use CracTools::SimCT::Fusion;
use CracTools::SimCT::Fusion::FusedExon;
use CracTools::SimCT::Annotations;
use CracTools::SimCT::GenomeSimulator;
use CracTools::Utils;
use Data::Dumper;
use Inline::Files 0.68;
use File::Temp;

# Load the fasta file
my $chr1_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
while(<CHR1_FASTA>) {print $chr1_fasta_file $_;}
close $chr1_fasta_file;
my $chr2_fasta_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
while(<CHR2_FASTA>) {print $chr2_fasta_file $_;}
close $chr2_fasta_file;

my $genome = CracTools::SimCT::Genome->new(
  reference_sequence_files => {
    1 => $chr1_fasta_file->filename,
    2 => $chr2_fasta_file->filename,
  },
);

# Load the GTF annotations
my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
while(<ANNOTATIONS>) {print $gtf_file $_;}
close $gtf_file;

my $annotations = CracTools::SimCT::Annotations->new();
$annotations->loadGTF($gtf_file->filename);

# Create the GenomeSimulator
my $gs = CracTools::SimCT::GenomeSimulator->new(
  genome => $genome,
);

$gs->addFusion(
  CracTools::SimCT::Fusion->new(
    fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(
      exon => $annotations->getGene('geneA')->getExon('11,27'),
    ),
    fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(
      exon => $annotations->getGene('geneB')->getExon('16,22'),
    ),
  ),
);

$gs->addFusion(
  CracTools::SimCT::Fusion->new(
    fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(
      exon => $annotations->getGene('geneB')->getExon('16,22'),
    ),
    fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(
      exon => $annotations->getGene('geneA')->getExon('11,27'),
    ),
  ),
);

$gs->addFusion(
  CracTools::SimCT::Fusion->new(
    fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(
      exon => $annotations->getGene('geneA')->getExon('11,27'),
    ),
    fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(
      exon => $annotations->getGene('geneA')->getExon('8,9'),
    ),
  ),
);

my $genome_dir = File::Temp->newdir();
my $gtf_output = new File::Temp( SUFFIX => '.gtf' );

my $sg = $gs->generateGenome(
  genome_dir  => $genome_dir,
  annotations => $annotations,
);

my $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/$CracTools::SimCT::Const::CHR_FUSIONS.fa");
my $entry     = $fasta_it->();


# CTAGCTAGTTAGCTCGATCGGC => GCCGATCGAGCTAACTAGCTAG
# lenght: 49
my $fusion_1 = "(TGGTAGTACCCGTCGCATGTCGAAAGT)(GCCGATCGAGCTAACTAGCTAG)";
# GATCGGC => GCCGATC
# length: 24
my $fusion_2 = "(GCCGATC)(CGTCGCATGTCGAAAGT)";
# length: 47
my $fusion_3 = "(TGGTAGTACCCGTCGCATGTCGAAAGT)(ACCCGTCGCATGTCGAAAGT)";
ok($entry->{seq} =~ /^($fusion_1)($fusion_2)($fusion_3)$/,"generateGenome - FASTA control (2)");

my $gtf_it        = CracTools::Utils::gffFileIterator("$genome_dir/annotations.gtf",'gtf');

# Skip regular exons
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

  my @alignments      = $sg->liftover->getAlignments($CracTools::SimCT::Const::CHR_FUSIONS,11,30);
  #print STDERR Dumper(\@alignments);
}

# Check the second fusion
{
  my @alignments      = $sg->liftover->getAlignments($CracTools::SimCT::Const::CHR_FUSIONS,50,70);
  #print STDERR Dumper(\@alignments);
}


# Check the third fusion
{
  my @alignments      = $sg->liftover->getAlignments($CracTools::SimCT::Const::CHR_FUSIONS,75,115);
  #print STDERR Dumper(\@alignments);
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

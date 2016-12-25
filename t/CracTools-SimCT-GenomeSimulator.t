use strict;
use warnings;

use Test::More tests => 54;
use CracTools::SimCT::Const;
use CracTools::SimCT::Utils;
use CracTools::SimCT::Genome;
use CracTools::SimCT::Annotations;
use CracTools::SimCT::GenomeSimulator;
use CracTools::SimCT::Fusion;
use CracTools::SimCT::Fusion::FusedExon;
use CracTools::SimCT::Mutation;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;
use CracTools::SimCT::GenomicInterval;
use CracTools::Utils;
use Inline::Files 0.68;
use File::Temp;

sub newInterval {
  my ($chr,$start,$end,$strand) = @_;
  return CracTools::SimCT::GenomicInterval->new(
    chr   => $chr,
    start => $start,
    end   => $end,
    strand => defined $strand? $strand : '+',
  );
}

{
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

  is($gs->allMutations,0);

  # Add mutations
  $gs->addMutation(
    CracTools::SimCT::Mutation::Substitution->new(
      chr => "1",
      start => 10,
      new_nuc => 'G',
    ),
  );
  $gs->addMutation(
    CracTools::SimCT::Mutation::Deletion->new(
      chr => "1",
      start => 15,
      length => 10,
    ),
  );
  $gs->addMutation(
    CracTools::SimCT::Mutation::Insertion->new(
      chr => "1",
      start => 2,
      inserted_sequence => 'AGG',
    ),
  );



  # chr1
  #                111111111----------1222222222233333333334444444444555
  #      0123456789012345678----------9012345678901234567890123456789012
  # sim: ATAGGGGTAGTACGCGTCA----------GTCT
  # ref: AT---GGTAGTACCCGTCGCATGTCGAAAGTCT
  #                   [ E1 ]    [ E2 ]
  #      01---234567890123456789012345678901234567890123456789012
  #                   1111111111222222222233333333334444444444555
  #
  # chr2
  #

  # Create a fusion gene with the first exon of geneA and the second exon of geneB
  my $gene_A  = $annotations->getGene('geneA');
  my @exons_A = $gene_A->sortedExons;
  my $exon_A_first_exon = shift @exons_A;
  my $exon_A_last_exon  = pop @exons_A;
  my $gene_B = $annotations->getGene('geneB');
  my @exons_B = $gene_B->sortedExons;
  my $exon_B_first_exon = shift @exons_B;
  my $exon_B_last_exon  = pop @exons_B;

  $gs->addFusion(CracTools::SimCT::Fusion->new(
      fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_A_first_exon),
      fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_B_last_exon),
  ));
  $gs->addFusion(CracTools::SimCT::Fusion->new(
      fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_A_last_exon),
      fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_A_first_exon),
  ));
  $gs->addFusion(CracTools::SimCT::Fusion->new(
      fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_B_first_exon),
      fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_A_first_exon),
  ));
  # Add a false fusion
  $gs->addFusion(CracTools::SimCT::Fusion->new(
      fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_A_first_exon),
      fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_A_last_exon),
  ));
  # Add a class 3 fusion
  $gs->addFusion(CracTools::SimCT::Fusion->new(
      fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_A_last_exon),
      fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_A_first_exon),
  ));

  my $genome_dir = File::Temp->newdir();
  #print STDERR "$genome_dir\n";

  # TODO sg should have a "genome" attribute
  my $sg = $gs->generateGenome(
    genome_dir   => $genome_dir,
    annotations  => $annotations,
  );
  #sleep 20;

  # Verify if fasta is good
  my $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/1.fa");
  my $entry     = $fasta_it->();
  is($entry->{seq},"ATAGGGGTAGTACGCGTCAGTCT","generateGenome - FASTA control (1)");

  {
    $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/fusion_0.fa");
    $entry     = $fasta_it->();
    # GeneA, exon1: CCCGTC
    # GeneB, exon1-exon2: GCTAGTTAGCTCGATC => GATCGAGCTAACTAGC
    my $fusion_1 = "CCCGTCGATCGAGCTAACTAGC";
    ok($entry->{seq} =~ /^$fusion_1$/,"generateGenome - FASTA control (2)");
  }
  {
    $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/fusion_1.fa");
    $entry     = $fasta_it->();
    # GeneA, exon1: CCCGTC
    # GeneB, exon1-exon2: GCTAGTTAGCTCGATC => GATCGAGCTAACTAGC
    my $fusion_2 = "(CCCGTCGCATGTCGAA){2}";
    ok($entry->{seq} =~ /^$fusion_2$/,"generateGenome - FASTA control (3)");
  }
  {
    $fasta_it  = CracTools::Utils::seqFileIterator("$genome_dir/fusion_2.fa");
    $entry     = $fasta_it->();
    # Check second fusion
    my $fusion_3 = "GATCGAGCTAACTAGCCCCGTCGCATGTCGAA";
    ok($entry->{seq} =~ /^$fusion_3$/,"generateGenome - FASTA control (4)");
  }

  # Verify liftover functions
  my $shifted_interval = $sg->liftover->shiftInterval(newInterval("1",15,20));
  is($shifted_interval->[0]->start,12,"shiftInterval (1)");
  is($shifted_interval->[0]->end,14,"shiftInterval (2)");
  is($shifted_interval->[1]->start,25,"shiftInterval (3)");
  is($shifted_interval->[1]->end,27,"shiftInterval (4)");

  # Verify liftover functions for fusions
  $shifted_interval = $sg->liftover->shiftInterval(newInterval("fusion_0",2,10));
  is($shifted_interval->[0]->chr,1,"shiftInterval (1)");
  is($shifted_interval->[0]->start,11,"shiftInterval (1)");
  is($shifted_interval->[0]->end,14,"shiftInterval (2)");
  is($shifted_interval->[0]->strand,'+',"shiftInterval (2)");
  is($shifted_interval->[1]->chr,2,"shiftInterval (1)");
  is($shifted_interval->[1]->start,15,"shiftInterval (1)");
  is($shifted_interval->[1]->end,19,"shiftInterval (2)");
  is($shifted_interval->[1]->strand,'-',"shiftInterval (2)");

  {
    my @alignments = $sg->liftover->getSplicedAlignments(
      newInterval("fusion_0",2,9),
      newInterval("fusion_0",12,17),
    );
    is($alignments[0]->chr,1,"chimeric spliced alignment");
    is($alignments[0]->start,11,"chimeric spliced alignment");
    is($alignments[0]->cigar,'4M10S',"chimeric spliced alignment");
    is($alignments[1]->chr,2,"chimeric spliced alignment");
    is($alignments[1]->start,8,"chimeric spliced alignment");
    is($alignments[1]->cigar,'6M2N4M4S',"chimeric spliced alignment");
    is($alignments[1]->strand,'-',"chimeric spliced alignment");
  }
  {
    my @alignments = $sg->liftover->getSplicedAlignments(
      newInterval("fusion_0",2,9,'-'),
      newInterval("fusion_0",12,17,'-'),
    );
    #use Data::Dumper;
    #print STDERR Dumper(\@alignments);
    is($alignments[0]->chr,2,"chimeric spliced alignment");
    is($alignments[0]->start,8,"chimeric spliced alignment");
    is($alignments[0]->cigar,'6M2N4M4S',"chimeric spliced alignment");
    is($alignments[0]->strand,'+',"chimeric spliced alignment");
    is($alignments[1]->chr,1,"chimeric spliced alignment");
    is($alignments[1]->start,11,"chimeric spliced alignment");
    is($alignments[1]->cigar,'10S4M',"chimeric spliced alignment");
    is($alignments[1]->strand,'-',"chimeric spliced alignment");
  }
  {
    my @alignments = $sg->liftover->getSplicedAlignments(
      newInterval("fusion_0",2,10,'+'),
    );
    is($alignments[0]->chr,1,"chimeric spliced alignment");
    is($alignments[0]->start,11,"chimeric spliced alignment");
    is($alignments[0]->cigar,'4M5S',"chimeric spliced alignment");
    is($alignments[1]->chr,2,"chimeric spliced alignment");
    is($alignments[1]->start,15,"chimeric spliced alignment");
    is($alignments[1]->cigar,'4S5M',"chimeric spliced alignment");
    is($alignments[1]->strand,'-',"chimeric spliced alignment");
  }
  {
    my @alignments = $sg->liftover->getSplicedAlignments(
      newInterval("fusion_0",2,10,'-'),
    );
    is($alignments[0]->chr,2,"chimeric spliced alignment");
    is($alignments[0]->start,15,"chimeric spliced alignment");
    is($alignments[0]->cigar,'5M4S',"chimeric spliced alignment");
    is($alignments[0]->strand,'+',"chimeric spliced alignment");
    is($alignments[1]->chr,1,"chimeric spliced alignment");
    is($alignments[1]->start,11,"chimeric spliced alignment");
    is($alignments[1]->cigar,'5S4M',"chimeric spliced alignment");
  }

  # Verify if annotations are good
  my $gtf_it        = CracTools::Utils::gffFileIterator("$genome_dir/annotations.gtf",'gtf');
  my $first_exon    = $gtf_it->();
  my $second_exon   = $gtf_it->();
  my $third_exon    = $gtf_it->();
  is($first_exon->{start},13,"generateGenome - GTF control (1)");
  is($first_exon->{end},18,"generateGenome - GTF control (2)");
  # Exon 2 is supressed by the deletion
  is($second_exon->{chr}, 2,"generateGenome - GTF control (3)");
  # Check the first fusion
  {
    my $fusion_exon_1   = $gtf_it->();
    my $fusion_exon_2   = $gtf_it->();
    is($fusion_exon_1->{chr},"fusion_0");
    is($fusion_exon_1->{start},1);
    is($fusion_exon_1->{end},12);
    is($fusion_exon_2->{start},17);
    is($fusion_exon_2->{end},22);
  }
  # TODO check VCF output
}

__CHR1_FASTA__
>chr1
ATGGTAGTACCCGTCGCATGTCGAA
AGTCT
__CHR2_FASTA__
>chr2
GCTAGCTAGTTAGCTCGATCGGCAT
__ANNOTATIONS__
1	StringTie	transcript	10	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	10	15	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	20	25	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	5	10	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";
2	StringTie	exon	15	20	1000	-	.	gene_id "geneB"; transcript_id "geneB.1"; exon_number "2"; cov "161.272552";

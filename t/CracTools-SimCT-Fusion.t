use strict;
use warnings;

use Test::More tests => 16;
use CracTools::SimCT::Fusion;
use CracTools::SimCT::Fusion::FusedExon;
use CracTools::SimCT::Annotations::Gene;
use CracTools::SimCT::Annotations::Exon;

#            [--------]                   [---]               [-------]
my $chr1 = "ATAGTCGACTAGCTATTAGATCGATCGGCGCTCTTTCGATCAGCTGAGTACGATCGACAT";
#                [----]             [-------]  
my $chr2 = "GCTGAGGCTTATATACCCAGTGTTCGGTAGCTAGTG";

# Gene A
my $gene_A = CracTools::SimCT::Annotations::Gene->new(
  chr     => 'chr1',
  strand  => '+',
  id      => 'gene_A',
);

my $exon_A1 = CracTools::SimCT::Annotations::Exon->new(
  start => 1,
  end   => 10,
  gene  => $gene_A,
);
my $exon_A2 = CracTools::SimCT::Annotations::Exon->new(
  start => 30,
  end   => 34,
  gene  => $gene_A,
);
my $exon_A3 = CracTools::SimCT::Annotations::Exon->new(
  start => 50,
  end   => 58,
  gene  => $gene_A,
);

# Gene B
my $gene_B = CracTools::SimCT::Annotations::Gene->new(
  chr     => 'chr2',
  strand  => '-',
  id      => 'gene_B',
);

my $exon_B1 = CracTools::SimCT::Annotations::Exon->new(
  start => 5,
  end   => 10,
  gene  => $gene_B,
);
my $exon_B2 = CracTools::SimCT::Annotations::Exon->new(
  start => 24,
  end   => 32,
  gene  => $gene_B,
);

{
  my $fused_exon_5 = CracTools::SimCT::Fusion::FusedExon::5prim->new(exon => $exon_A2);
  my $fused_exon_3 = CracTools::SimCT::Fusion::FusedExon::3prim->new(exon => $exon_B2);

  $fused_exon_5->setFusedSequence(\$chr1);
  $fused_exon_3->setFusedSequence(\$chr2);

  is($fused_exon_5->fused_sequence,"TAGTCGACTAGCTATTAGATCGATCGGCGCTCTT");
  is($fused_exon_3->fused_sequence,"TAGCTACCGAACACTGGGTATATAAGCC");
  
  my @fused_exons_5 = $fused_exon_5->allFusedExons;
  is(scalar @fused_exons_5,2);
  my $exon_test = pop @fused_exons_5;
  is($exon_test->start,30);
  is($exon_test->end,34);

  my @fused_exons_3 = $fused_exon_3->allFusedExons;
  is(scalar @fused_exons_3,2);
  $exon_test = shift @fused_exons_3;
  is($exon_test->start,24);
  is($exon_test->end,32);

  my $fusion = CracTools::SimCT::Fusion->new(
    fused_exon_5prim  => $fused_exon_5,
    fused_exon_3prim  => $fused_exon_3,
  );

  $fusion->chr_fusion_pos(0);

  my $fusion_gene = $fusion->getFusionGene;
  is($fusion_gene->start,0);
  is($fusion_gene->end,61);

  my @fusion_exons = $fusion_gene->sortedExons;
  is(scalar @fusion_exons, 3);
  $exon_test = shift @fusion_exons;
  is($exon_test->end,9);
  $exon_test = shift @fusion_exons;
  is($exon_test->start,29);
  is($exon_test->end,42);
  $exon_test = shift @fusion_exons;
  is($exon_test->start,56);
  is($exon_test->end,61);
}

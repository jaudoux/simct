package CracTools::SimCT::MutationGenerator::Fusion;  
# ABSTRACT: A mutation generator that introduce random fusions

use Moose;

use CracTools::Const;
use CracTools::SimCT::Const;
use CracTools::SimCT::Fusion;
use CracTools::SimCT::Fusion::FusedExon;

#has annotations => (
#  is => 'ro',
#  isa => 'CracTools::SimCT::Annotations',
#  required => 1,
#);

with 'CracTools::SimCT::MutationGenerator';

# Generate random mutations in the genomeSimulator object
sub generateMutations {
  my $self        = shift;
  my $annotations = shift;
  my $nb_fusions  = shift;

  # If there is not 2 genes at least, then we can not create
  # any fusions
  return 0 if $annotations->nbGenes < 2;

  while($nb_fusions > 0) {
    # We pick to random genes out of the gene catalogue
    my $gene_A = ($annotations->allGenes)[int rand $annotations->nbGenes];
    my $gene_B = ($annotations->allGenes)[int rand $annotations->nbGenes];

    # If we have pick the same gene two times we try again
    #next if $gene_A eq $gene_B;

    # If the fusion is properly added, we decrement our counter
    $nb_fusions-- if $self->genome_simulator->addFusion(
      CracTools::SimCT::Fusion->new(
        fused_exon_5prim => CracTools::SimCT::Fusion::FusedExon::5prim->new(
          exon => ($gene_A->allExons)[int rand $gene_A->nbExons],
        ),
        fused_exon_3prim => CracTools::SimCT::Fusion::FusedExon::3prim->new(
          exon => ($gene_B->allExons)[int rand $gene_B->nbExons],
        ),
      ),
    );
  }

  return 1;
}

1;

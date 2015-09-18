use strict;
use warnings;

package CracTools::SimCT::MutationGenerator::Fusion;  
# ABSTRACT: A mutation generator that introduce random fusions

use parent 'CracTools::SimCT::MutationGenerator';

use Carp;

use CracTools::Const;
use CracTools::SimCT::Const;

sub _init {
  my $self = shift;
  my $args = shift;
  # TODO Add a "min_fusion_distance" parameter!!!
}

# Generate random mutations in the genomeSimulator object
sub generateMutations {
  my $self        = shift;
  my $nb_fusions  = shift;

  # If there is not 2 genes at least, then we can not create
  # any fusions
  return 0 if $self->genomeSimulator->genes < 2;

  while($nb_fusions > 0) {
    # We pick to random genes out of the gene catalogue
    my $gene_A = ($self->genomeSimulator->genes)[int rand $self->genomeSimulator->genes];
    my $gene_B = ($self->genomeSimulator->genes)[int rand $self->genomeSimulator->genes];

    # If we have pick the same gene two times we try again
    next if $gene_A eq $gene_B;

    # If the fusion is properly added, we decrement our counter
    $nb_fusions-- if $self->genomeSimulator->addFusion(
      $gene_A,
      # We choose one random exon of gene_A
      ($self->genomeSimulator->exons($gene_A))[int rand $self->genomeSimulator->exons($gene_A)],
      $gene_B,
      # We choose one random exon of gene_B
      ($self->genomeSimulator->exons($gene_B))[int rand $self->genomeSimulator->exons($gene_B)],
    );
  }

  return 1;
}

1;

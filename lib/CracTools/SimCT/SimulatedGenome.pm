package CracTools::SimCT::SimulatedGenome;
# ABSTRACT: Represent a Simulated genome and contains mapping values between the original genome

use Moose;

require 5.004;
use Tie::RefHash;

use CracTools::Utils;
use CracTools::SimCT::Utils;
use CracTools::SimCT::LiftOver;
use CracTools::SimCT::MutationQuery;

has 'liftover' => (
  is  => 'ro',
  isa => 'CracTools::SimCT::LiftOver',
  default => sub { CracTools::SimCT::LiftOver->new(); },
);

has 'mutation_query' => (
  is  => 'ro',
  isa => 'CracTools::SimCT::MutationQuery',
  default => sub { CracTools::SimCT::MutationQuery->new(); },
);

sub BUILD {
  my $self = shift;
  my $args = shift;
  my $genome_simulator = $args->{genome_simulator};

  # Create a liftover object
  foreach my $chr (sort $genome_simulator->genome->references) {
    my $prev_pos        = 0;
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome
    my $chr_length      = $genome_simulator->genome->getReferenceLength($chr);

    # FIXME the liftover build here, should only support sub / ins / del, but not
    # complexe mutations. Look into GenomeSimulator for a valid examples.
    foreach my $mut ($genome_simulator->sortedMutations($chr)) {

      # Get the mutations pos on the simulated genome
      my $mut_sg_pos = $mut->start - $offset;

      # add the mutation into the mutation query
      $self->mutation_query->addMutation($mut);

      # Skip mutations that have no incidence on the offset
      next if $mut->referenceLength == 1 && $mut->mutationLength == 1;

      # Update prev_pos and index
      $prev_pos = $index;
      $index    = $mut_sg_pos;

      # Add the current interval
      $self->liftover->addInterval($chr,$prev_pos,$index-1,$offset);

      # Update offset and index
      $index  += $mut->mutationLength;
      $offset += ($mut->referenceLength - $mut->mutationLength);
    }
    # Add the last interval
    if($index < $chr_length) {
      $self->liftover->addInterval($chr,$index,$chr_length-1,$offset);
    }
  }

  # Add fusions in the liftover query
  foreach my $fusion ($genome_simulator->allFusions) {
    my $start_pos   = $fusion->chr_fusion_pos;
    foreach my $fused_exon (($fusion->fused_exon_5prim,$fusion->fused_exon_3prim)) {
      # Add interval for the 5prim part
      $self->liftover->addInterval(
        $fusion->chr_fusion,
        $start_pos,
        $start_pos + $fused_exon->length - 1,
        $fused_exon->getFusionOffset($start_pos),
        $fused_exon->chr,
        $fused_exon->strand eq '+'? 0 : 1,
      );
      # Update the start pos for 3prim fused exon
      $start_pos    += $fused_exon->length;
    }
  }
};

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 DESCRIPTION


=head1 ACCESSORS

=head2 liftover

=head2 mutation_query

=head1 METHODS

=head2 getAlignements

    Arg[1]  :  'CracTools::SimCT::GenomicInterval' - genomic interval

Given a genomic interval over the simulated genome, returns a list of tha alignement
over the reference genome as a list of genomic intervals.

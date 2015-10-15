package CracTools::SimCT::SimulatedGenome;
# ABSTRACT: Represent a Simulated genome and contains mapping values between the original genome

use Moose;

require 5.004;
use Tie::RefHash;

use CracTools::Utils;
use CracTools::SimCT::Utils;

has 'liftover' => (
  is  => 'ro',
  isa => 'CracTools::SimCT::LiftOver',
  default => sub { CracTools::SimCT::LiftOver->new(); },
);

has 'mutation_query' => (
  is  => 'ro',
  isa => 'CracTools::Interval::Query',
  default => sub { CracTools::Interval::Query->new(); },
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

    foreach my $mut ($genome_simulator->sortedMutations($chr)) {
      
      #next if ($mut->referenceLength - $mut->mutationLength) == 0;

      $prev_pos = $index;
      $index    = $mut->pos - $offset;

      # Add the current interval
      $self->liftover->addInterval($chr,$prev_pos,$index-1,$offset);

      # Add the mutation into the mutation query
      $self->mutation_query->addInterval($chr,$index,$index,1,$mut);

      # Update offset and index
      $index  += $mut->mutationLength; 
      $offset += ($mut->referenceLength - $mut->mutationLength);
    }
    # Add the last interval
    $self->liftover->addInterval($chr,$index,$chr_length-1,$offset);
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

1;

__END__

=head1 DESCRIPTION


=head1 ACCESSORS

=head2 liftover

=head2 mutation_query

=head1 METHODS

package CracTools::SimCT::LiftOver::ShiftedInterval;
# ABSTRACT: Represents a chunk of a genomic intervals that have been lifted over a different refernce

use Moose;

extends 'CracTools::SimCT::GenomicInterval';

use CracTools::SimCT::GenomicInterval;

has 'reference_interval' => (
  isa => 'CracTools::SimCT::GenomicInterval',
  is  => 'ro',
);

sub isReversed {
  my $self = shift;
  return $self->strand ne $self->reference_interval->strand;
}

no Moose;
__PACKAGE__->meta->make_immutable;

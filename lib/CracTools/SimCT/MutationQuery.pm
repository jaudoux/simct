package CracTools::SimCT::MutationQuery;
# ABSTRACT: A mutation query object holds a set of mutations and provide a way of querying them

use Moose;

use CracTools::Utils;

extends 'CracTools::Interval::Query';

sub addMutation {
  my $self  = shift;
  my $mut   = shift;

  # add the mutation into the mutation query
  $self->addInterval($mut->chr,$mut->start,$mut->end,1,$mut);
}

sub addMutations {
  my $self  = shift;
  my $mutations = shift;

  foreach my $mut (@{$mutations}) {
    $self->addMutation($mut);
  }
}

sub getOverlappingMutations {
  my $self = shift;
  my @intervals = @_;
  my @mutations;

  foreach my $interval (@intervals) {
    push @mutations, @{$self->fetchByRegion(
      $interval->chr,
      $interval->start,
      $interval->end,
    )};
  }
  return @mutations;
}

no Moose;
# Not inlining the constructor because CracTools::Interval::Query is not moose-based
__PACKAGE__->meta->make_immutable(inline_constructor => 0);

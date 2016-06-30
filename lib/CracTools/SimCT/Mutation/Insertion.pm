package CracTools::SimCT::Mutation::Insertion;
# ABSTRACT: Mutation of type 'insertion'

use Moose;

extends 'CracTools::SimCT::Mutation';

sub referenceLength {
  my $self = shift;
  return 0;
}

sub mutation_sequence {
  my $self = shift;
  return $self->inserted_sequence;
}

#with 'CracTools::SimCT::Mutation';

has inserted_sequence => (
  is  => 'ro',
  isa => 'DNA',
  required => 1,
);

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 METHODS

=head2 new

  Arg [chr]               : 'Str'   - insertion's chromosome
  Arg [start]               : 'Int'   - insertion's strand
  Arg [inserted_sequence] : 'Str'   - inserted sequence

Create a new 'CracTools::SimCT::Mutation::Insertion'

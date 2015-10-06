package CracTools::SimCT::Mutation::Insertion;

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


1;

__END__

=head1 METHODS

=head2 new

  Arg [chr]               : 'Str'   - insertion's chromosome
  Arg [pos]               : 'Int'   - insertion's strand
  Arg [inserted_sequence] : 'Str'   - inserted sequence

Create a new 'CracTools::SimCT::Mutation::Insertion'

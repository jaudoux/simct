package CracTools::SimCT::Mutation::Deletion;
# ABSTRACT: Mutation of type 'deletion'

use Moose;

extends 'CracTools::SimCT::Mutation';

sub referenceLength {
  my $self = shift;
  return $self->length;
}

sub mutation_sequence {
  my $self = shift;
  return "";
}

#with 'CracTools::SimCT::Mutation';

has length => (
  is        => 'ro',
  isa       => 'Natural',
  required  => 1,
);

1;

__END__

=head1 METHODS

=head2 new

  Arg [chr]    : 'Str'   - deletion's chromosome
  Arg [pos]    : 'Int'   - deletion's strand
  Arg [length] : 'Int'   - deletion length

Create a new 'CracTools::SimCT::Mutation::Deletion'

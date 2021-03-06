package CracTools::SimCT::Mutation::Substitution;
# ABSTRACT: Mutation of type 'substitution'

use Moose;

extends 'CracTools::SimCT::Mutation';

sub referenceLength {
  my $self = shift;
  return 1;
}

sub mutation_sequence {
  my $self = shift;
  return $self->new_nuc;
}

#with 'CracTools::SimCT::Mutation';

has new_nuc => (
  is  => 'ro',
  isa => 'DNAnuc',
  required => 1,
);

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 METHODS

=head2 new

  Arg [chr]     : 'Str'   - substitution's chromosome
  Arg [start]     : 'Int'   - substitution's strand
  Arg [new_nuc] : 'Str'   - substitution nucleotide

Create a new 'CracTools::SimCT::Mutation::Substitution'

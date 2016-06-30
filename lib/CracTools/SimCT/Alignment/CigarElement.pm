package CracTools::SimCT::Alignment::CigarElement;
# ABSTRACT: A cigar element to describe an alignement

use Moose;

use CracTools::SimCT::Utils;

has 'op'  => (
  isa => 'CigarOperator',
  is  => 'rw',
);

has 'nb'  => (
  isa => 'Natural',
  is  => 'rw',
);

sub isReferenceBased {
  my $self = shift;
  return $self->op =~ /^[M=XDN]$/;
}

no Moose;
__PACKAGE__->meta->make_immutable;

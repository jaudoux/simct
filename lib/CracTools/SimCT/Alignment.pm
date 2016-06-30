package CracTools::SimCT::Alignment;
# ABSTRACT: An alignment

use Moose;

extends 'CracTools::SimCT::GenomicInterval';

use CracTools::SimCT::Alignment::CigarElement;

has 'cigar' => (
  traits  => ['Array'],
  is      => 'ro',
  isa     => 'ArrayRef[CracTools::SimCT::Alignment::CigarElement]',
  default => sub { [] },
  reader => '_cigar',
  handles => {
    allCigarElements    => 'elements',
    appendCigarElement  => 'push',
    prependCigarElement => 'unshift',
  },
);

has '+end'  => (
  lazy  => 1,
  init_arg => undef,
  default => sub {
    my $self = shift;
    return $self->start;
  },
);

around 'end' => sub {
  my $orig = shift;
  my $self = shift;
  return $self->start + $self->_sumCigel('M','D','N','X','=') - 1;
};

has 'query_mapping_start' => (
  is      => 'rw',
  isa     => 'Int',
  default => 0,
);

has 'query_length' => (
  is        => 'rw',
  isa       => 'Int',
  # It could not be defined, and computed from the query mapping length
  required  => 1,
  # We should do something to verify that the query length is never lower
  # than the query mapping length.
  #trigger => \&_query_length_set,
);

sub query_mapping_length {
  my $self = shift;
  return $self->_sumCigel('M','X','=','I');
}

sub query_mapping_end {
  my $self = shift;
  return $self->query_mapping_start + $self->query_mapping_length - 1;
}

sub right_softclip {
  my $self = shift;
  return $self->query_length - ($self->query_mapping_start + $self->query_mapping_length);
}

sub cigar {
  my $self = shift;
  my $cigar = join "", map { $_->nb.$_->op } $self->allCigarElements;

  # append softclip if needed
  #my $left_softclip = $self->query_mapping_start;
  if($self->query_mapping_start > 0) {
    $cigar = $self->query_mapping_start."S".$cigar;
  }

  # prepend softclip if needed
  if($self->right_softclip > 0) {
    $cigar .= $self->right_softclip."S";
  }
  return $cigar;
}

sub _sumCigel {
  my $self = shift;
  my @cigar_operators = @_;
  my $sum = 0;
  my $regex = "(";
  $regex .= join "|", @cigar_operators;
  $regex .= ")";
  foreach my $cigel ($self->allCigarElements) {
    if($cigel->op =~ /$regex/) {
      $sum += $cigel->nb;
    }
  }
  return $sum;
}

no Moose;
__PACKAGE__->meta->make_immutable;

package CracTools::SimCT::Mutation;
# ABSTRACT: Base class for mutations

#use Moose::Role;
use Moose;

extends 'CracTools::SimCT::GenomicInterval';

use CracTools::SimCT::Utils;

has '+end' => (
  lazy    => 1,
  default => sub {
    my $self = shift;
    if($self->referenceLength == 0) {
      return $self->start;
    } else {
      return $self->start + $self->referenceLength - 1;
    }
  },
);

has 'reference_sequence' => (
  is        => 'rw',
  isa       => 'DNA',
  # By default we set the reference sequence to N's
  lazy      => 1,
  default   => sub {
    my $self = shift;
    my $seq  = "";
    $seq .= 'N' for(1..$self->referenceLength);
    return $seq;
  },
);

has 'frequency' => (
  is        => 'rw',
  isa       => 'Num',
  default   => 1,
);

sub mutation_sequence {
  #my $self = shift;
  return "";
}

sub referenceLength {
  return 0;
}

sub mutationLength {
  my $self = shift;
  return length $self->mutation_sequence;
}

sub setReferenceSequence {
  my $self = shift;
  my $ref_seq = shift;
  die "Reference sequence has not the good length"
    unless length $ref_seq == length $self->reference_sequence;
  $self->reference_sequence($ref_seq);
}

sub getVCFRecord {
  my $self = shift;
  my %vcf_record = (
    chr => $self->chr,
    pos => $self->start,
    ref => $self->reference_sequence,
    alt => $self->mutation_sequence,
    id  => '.',
  );
  if($self->referenceLength == 0 || $self->mutationLength == 0) {
    $vcf_record{pos}--;
    $vcf_record{ref} = "N".$vcf_record{ref};
    $vcf_record{alt} = "N".$vcf_record{alt};
  }
  return \%vcf_record;
}

#with ('CracTools::SimCT::GeneticVariant');

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 SYNOPSIS

=head1 ATTRIBUTES

=head2 chr

The chromosome where the mutation occurs

=head2 start

The position of the mutation on the reference genome

=head2 reference_sequence

The sequence of the mutation on the reference genome.

=head2 mutation_sequence

The sequence of the mutation on the mutated genome

=head1 METHODS

=head2 referenceLength

The length of the mutation on the reference genome

=head2 mutationLength

The length of the mutation on the mutated genome

=head2 getVCFRecord => HashRef[chr => '', pos => '', alt => '', ref => '']

Return a hash reference of the VCF record generated for this mutation

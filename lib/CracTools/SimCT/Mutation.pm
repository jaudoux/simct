package CracTools::SimCT::Mutation;

#use Moose::Role;
use Moose;
use CracTools::SimCT::Utils;

has 'chr' => (
  is        => 'ro',
  isa       => 'Str',
  required  => 1,
);

has 'pos' => (
  is        => 'ro',
  isa       => 'Int',
  required  => 1,
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

sub mutation_sequence {
  #my $self = shift;
  return "";
}

sub referenceLength {
  return 0;
}

#requires 'referenceLength';

#requires 'mutation_sequence';

sub end {
  my $self = shift;
  if($self->referenceLength == 0) {
    return $self->pos;
  } else {
    return $self->pos + $self->referenceLength - 1;
  }
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
    pos => $self->pos,
    ref => $self->reference_sequence,
    alt => $self->mutation_sequence,
  );
  if($self->referenceLength == 0 || $self->mutationLength == 0) {
    $vcf_record{pos}--;
    $vcf_record{ref} = "N".$vcf_record{ref};
    $vcf_record{pos} = "N".$vcf_record{pos};
  }
  return \%vcf_record;
}

#with ('CracTools::SimCT::GeneticVariant');

1;

__END__

=head1 SYNOPSIS

=head1 ATTRIBUTES

=head2 chr

The chromosome where the mutation occurs

=head2 pos

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

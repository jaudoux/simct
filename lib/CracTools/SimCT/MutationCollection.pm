package CracTools::SimCT::MutationCollection;
# ABSTRACT: A mutation collection object holds a set of mutations and is able to load or write them.

use Moose;

use CracTools::Utils;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;

has 'substitutions' => (
  is  => 'ro',
  isa => 'ArrayRef[CracTools::SimCT::Mutation::Substitution]',
  default => sub { [] },
);

has 'deletions' => (
  is  => 'ro',
  isa => 'ArrayRef[CracTools::SimCT::Mutation::Deletion]',
  default => sub { [] },
);

has 'insertions' => (
  is  => 'ro',
  isa => 'ArrayRef[CracTools::SimCT::Mutation::Insertion]',
  default => sub { [] },
);

sub loadVCF {
  my $self      = shift;
  my $vcf_file  = shift;
  my $genome    = shift;
  my $max_ins   = shift;
  my $max_del   = shift;

  my $vcf_it = CracTools::Utils::vcfFileIterator($vcf_file);
  while(my $vcf_record = $vcf_it->()) {
    my $chr = $vcf_record->{chr};
    # Do not load the mutations if the corresponding reference
    # is not available in the genome
    next if defined $genome && !defined $genome->getReferenceFile($chr);

    my $ref_length = length $vcf_record->{ref};
    my @AFs = defined $vcf_record->{info}->{AF}?
      split ',', $vcf_record->{info}->{AF} : ();

    # Loop over alternative and add a new mutation for each
    foreach my $alt (@{$vcf_record->{alt}}) {
      # Skip monomorphic mutations
      next if $alt eq '.';

      # Get alt informations
      my $alt_length = length $alt;
      my $frequency = shift @AFs;

      # SNP case
      if($ref_length == $alt_length) {
        my $sub_index = $ref_length - 1;
        push @{$self->substitutions}, CracTools::SimCT::Mutation::Substitution->new(
          chr => $chr,
          start => $vcf_record->{pos} + $sub_index,
          new_nuc => substr($alt, $sub_index, 1),
          frequency => defined $frequency? $frequency : 1,
        );

      # Insertion case
      } elsif($alt_length > $ref_length) {
        my $ins_length = $alt_length - $ref_length;
        # Skip insertions larger that the "max_ins"
        next if defined $max_ins && $ins_length > $max_ins;
        my $ins_index  = $alt_length - $ins_length;
        push @{$self->insertions}, CracTools::SimCT::Mutation::Insertion->new(
          chr => $chr,
          start => $vcf_record->{pos} + $ins_index,
          inserted_sequence => substr($alt, $ins_index),
          frequency => defined $frequency? $frequency: 1,
        );

      # Deletion case
      } else {
        my $del_length = $ref_length - $alt_length;
        # Skip deletions larger that the "max_del"
        next if defined $max_del && $del_length > $max_del;
        my $del_index  = $ref_length - $del_length;
        push @{$self->deletions}, CracTools::SimCT::Mutation::Deletion->new(
          chr => $chr,
          start => $vcf_record->{pos} + $del_index,
          length => $del_length,
          frequency => defined $frequency? $frequency : 1,
        );
      }
    }
  }
}

no Moose;
__PACKAGE__->meta->make_immutable;

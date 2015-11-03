package CracTools::SimCT::MutationGenerator::VCF;
# ABSTRACT: A mutation generator that introduce mutations from a VCF source

use Moose;

use List::Util qw(shuffle min);

use CracTools::Utils;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;

extends 'CracTools::SimCT::MutationGenerator::Random';

has 'substitutions' => (
  is  => 'ro',
  isa => 'HashRef[ArrayRef[CracTools::SimCT::Mutation::Substitution]]',
  default => sub { {} },
);

has 'deletions' => (
  is  => 'ro',
  isa => 'HashRef[ArrayRef[CracTools::SimCT::Mutation::Deletion]]',
  default => sub { {} },
);

has 'insertions' => (
  is  => 'ro',
  isa => 'HashRef[ArrayRef[CracTools::SimCT::Mutation::Insertion]]',
  default => sub { {} },
);

sub loadVCF {
  my $self = shift;
  my $vcf_file = shift;
  
  my $vcf_it = CracTools::Utils::vcfFileIterator($vcf_file);
  while(my $vcf_record = $vcf_it->()) {
    my $chr = $vcf_record->{chr};
    # Do not load the mutations if the corresponding reference
    # is not available in the genome
    next if !defined $self->genome_simulator->genome->getReferenceFile($chr);
    my $ref_length = length $vcf_record->{ref};
    foreach my $alt (@{$vcf_record->{alt}}) {
      my $alt_length = length $alt;
      # SNP case
      if($ref_length == $alt_length) {
        my $sub_index = $ref_length - 1;
        push @{$self->substitutions->{$chr}}, CracTools::SimCT::Mutation::Substitution->new(
          chr => $chr,
          pos => $vcf_record->{pos} + $sub_index,
          new_nuc => substr($alt, $sub_index, 1),
        );
      # Insertion case
      } elsif($alt_length > $ref_length) {
        my $ins_length = $alt_length - $ref_length;
        # Skip insertions larger that the "max_ins"
        next if $ins_length > $self->max_ins;
        my $ins_index  = $alt_length - $ins_length;
        push @{$self->insertions->{$chr}}, CracTools::SimCT::Mutation::Insertion->new(
          chr => $chr,
          pos => $vcf_record->{pos} + $ins_index,
          inserted_sequence => substr($alt, $ins_index),
        );
      # Deletion case
      } else {
        my $del_length = $ref_length - $alt_length;
        # Skip deletions larger that the "max_del"
        next if $del_length > $self->max_del;
        my $del_index  = $ref_length - $del_length;
        push @{$self->deletions->{$chr}}, CracTools::SimCT::Mutation::Deletion->new(
          chr => $chr,
          pos => $vcf_record->{pos} + $del_index,
          length => $del_length,
        );
      }
    }
  }
}

sub generateMutations {
  my $self = shift;

  # Then we add a sample of the loaded mutations depending on the rates
  foreach my $chr ($self->genome_simulator->genome->references) {
    my $chr_length  = $self->genome_simulator->genome->getReferenceLength($chr);

    # shuffle the mutations for this chr
    my @snps       = defined $self->substitutions->{$chr}? 
                      shuffle @{$self->substitutions->{$chr}} : ();
    my @insertions = defined $self->insertions->{$chr}? 
                      shuffle @{$self->insertions->{$chr}}    : ();
    my @deletions  = defined $self->deletions->{$chr}? 
                      shuffle @{$self->deletions->{$chr}}     : ();

    my $nb_sub      = min(int($chr_length * $self->sub_rate / 100),scalar @snps);
    my $nb_ins      = min(int($chr_length * $self->ins_rate / 100),scalar @insertions);
    my $nb_del      = min(int($chr_length * $self->del_rate / 100),scalar @deletions);

    $self->genome_simulator->addMutation($snps[$_])       for 0..$nb_sub-1;
    $self->genome_simulator->addMutation($insertions[$_]) for 0..$nb_ins-1;
    $self->genome_simulator->addMutation($deletions[$_])  for 0..$nb_del-1;
  }

  return 1;
}

1;

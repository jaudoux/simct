package CracTools::SimCT::MutationGenerator::VCF;
# ABSTRACT: A mutation generator that introduce mutations from a VCF source

use Moose;

use List::Util qw(shuffle min);

use CracTools::Utils;
use CracTools::SimCT::MutationCollection;

extends 'CracTools::SimCT::MutationGenerator::Random';

has 'mutation_collection' => (
  is  => 'ro',
  isa =>  'CracTools::SimCT::MutationCollection',
  default => sub { CracTools::SimCT::MutationCollection->new() },
);

sub loadVCF() {
  my $self = shift;
  my $vcf_file = shift;
  $self->mutation_collection->loadVCF($vcf_file,
    $self->genome_simulator->genome,
    $self->max_ins,
    $self->max_del,
  );
}

sub generateMutations {
  my $self = shift;
  my $vcf_file = shift;

  # Then we add a sample of the loaded mutations depending on the rates
  my $genome_length  = $self->genome_simulator->genome->genomeLength;

  # shuffle the mutations for this chr
  my @snps       = defined $self->mutation_collection->substitutions?
                    shuffle @{$self->mutation_collection->substitutions} : ();
  my @insertions = defined $self->mutation_collection->insertions?
                    shuffle @{$self->mutation_collection->insertions}    : ();
  my @deletions  = defined $self->mutation_collection->deletions?
                    shuffle @{$self->mutation_collection->deletions}     : ();

  my $nb_sub      = min(int($genome_length * $self->sub_rate / 100),scalar @snps);
  my $nb_ins      = min(int($genome_length * $self->ins_rate / 100),scalar @insertions);
  my $nb_del      = min(int($genome_length * $self->del_rate / 100),scalar @deletions);

  $self->genome_simulator->addMutation($snps[$_])       for 0..$nb_sub-1;
  $self->genome_simulator->addMutation($insertions[$_]) for 0..$nb_ins-1;
  $self->genome_simulator->addMutation($deletions[$_])  for 0..$nb_del-1;

  return 1;
}

1;

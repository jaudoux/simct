use strict;
use warnings;

package CracTools::SimCT::MutationGenerator::Random;  
# ABSTRACT: A mutation generator that introduce random mutations ('ins','del','sub')

use parent 'CracTools::SimCT::MutationGenerator';

use Carp;

use CracTools::Const;
use CracTools::SimCT::Const;

sub _init {
  my $self = shift;
  my %args = @_;
  
  # Mutation rates
  $self->{mutation_rates} = {
    sub => $CracTools::SimCT::Const::SUB_RATE,
    ins => $CracTools::SimCT::Const::INS_RATE,
    del => $CracTools::SimCT::Const::DEL_RATE,
  };
  # If user define its mutation rate we update them
  $self->setMutationRate('ins',$args{ins_rate}) if defined $args{ins_rate};
  $self->setMutationRate('del',$args{del_rate}) if defined $args{del_rate};
  $self->setMutationRate('sub',$args{sub_rate}) if defined $args{sub_rate};

  # Insertion and deletion length
  $self->{max_ins} = defined $args{max_ins}? $args{max_ins} : $CracTools::SimCT::Const::MAX_INS;
  $self->{max_del} = defined $args{max_del}? $args{max_del} : $CracTools::SimCT::Const::MAX_DEL;
}

sub setMutationRate {
  my $self = shift;
  my ($mutation_type,$mutation_rate) = @_;
  if(defined $mutation_type && defined $mutation_rate) {
    $self->{mutation_rates}->{$mutation_type} = $mutation_rate;
  } else {
    carp "Missing argument to set mutation rate";
  }
}

sub getMutationRate {
  my $self          = shift;
  my $mutation_type = shift;
  return $self->{mutation_rates}->{$mutation_type};
}

sub maxInsertion {
  my $self = shift;
  return $self->{max_ins};
}

sub maxDeletion {
  my $self = shift;
  return $self->{max_del};
}

# Generate random mutations in the genomeSimulator object
sub generateMutations {
  my $self = shift;
  foreach my $chr ($self->genomeSimulator->references) {
    my $chr_length  = $self->genomeSimulator->getReferenceLength($chr);
    my $nb_sub      = int($chr_length * $self->getMutationRate('sub') / 100);
    my $nb_ins      = int($chr_length * $self->getMutationRate('ins') / 100);
    my $nb_del      = int($chr_length * $self->getMutationRate('del') / 100);

    # Generate Substitution
    while($nb_sub > 0) {
      $nb_sub-- if $self->genomeSimulator->addSubstitution(
        $chr,
        int(rand($chr_length)),
        $CracTools::Const::NUCLEOTIDES->[int(rand(4))],
      );
    }

    # Generate Insertions
    while($nb_ins > 0) {
      $nb_ins-- if $self->genomeSimulator->addInsertion(
        $chr,
        int(rand($chr_length)),
        join("", map{$CracTools::Const::NUCLEOTIDES->[int(rand(4))]} (1..int(rand(10-1)+1))), # This generates a random insertion
      );
    }

    # Generate Deletions
    while($nb_del > 0) {
      $nb_del-- if $self->genomeSimulator->addDeletion(
        $chr,
        int(rand($chr_length - $self->maxDeletion)),
        int(rand($self->maxDeletion-1)) + 1,
      );
    }
  }
}

1;

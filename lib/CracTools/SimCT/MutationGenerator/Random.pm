package CracTools::SimCT::MutationGenerator::Random;  
# ABSTRACT: A mutation generator that introduce random mutations ('ins','del','sub')

use Moose;

use CracTools::Const;
use CracTools::SimCT::Const;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;

# Mutation rates
has 'ins_rate' => (
  is => 'rw',
  isa => 'Num',
  default => $CracTools::SimCT::Const::INS_RATE,
);
has 'del_rate' => (
  is => 'rw',
  isa => 'Num',
  default => $CracTools::SimCT::Const::DEL_RATE,
);
has 'sub_rate' => (
  is => 'rw',
  isa => 'Num', 
  default => $CracTools::SimCT::Const::SUB_RATE,
);

# Maximum insertion/deletion length
has 'max_ins' => (
  is => 'rw',
  isa => 'Int',
  default => $CracTools::SimCT::Const::MAX_INS,
);
has 'max_del' => (
  is => 'rw',
  isa => 'Int',
  default => $CracTools::SimCT::Const::MAX_DEL,
);

with 'CracTools::SimCT::MutationGenerator';

# Generate random mutations in the genomeSimulator object
sub generateMutations {
  my $self = shift;
  foreach my $chr ($self->genome_simulator->genome->references) {
    my $chr_length  = $self->genome_simulator->genome->getReferenceLength($chr);

    my $nb_sub      = int($chr_length * $self->sub_rate / 100);
    my $nb_ins      = int($chr_length * $self->ins_rate / 100);
    my $nb_del      = int($chr_length * $self->del_rate / 100);

    # Generate Substitution
    while($nb_sub > 0) {
      $nb_sub-- if $self->genomeSimulator->addMutation(
        CracTools::SimCT::Mutation::Substitution->new(
          chr => $chr,
          pos => int(rand($chr_length)),
          new_nuc => $CracTools::Const::NUCLEOTIDES->[int(rand(4))],
        ),
      );
    }

    # Generate Insertions
    while($nb_ins > 0) {
      $nb_ins-- if $self->genomeSimulator->addInsertion(
        CracTools::SimCT::Mutation::Insertion->new(
          chr => $chr,
          pos => int(rand($chr_length)),
          inserted_sequence => join("", 
            map{$CracTools::Const::NUCLEOTIDES->[int(rand(4))]} (1..int(rand(10-1)+1))
          ), # This generates a random insertion
        ),
      );
    }

    # Generate Deletions
    while($nb_del > 0) {
      $nb_del-- if $self->genomeSimulator->addDeletion(
        CracTools::SimCT::Mutation::Deletion->new(
          chr => $chr,
          pos => int(rand($chr_length - $self->maxDeletion)),
          length => int(rand($self->maxDeletion-1)) + 1,
        ),
      );
    }
  }
}

1;

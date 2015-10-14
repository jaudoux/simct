package CracTools::SimCT::MutationGenerator::Random;  
# ABSTRACT: A mutation generator that introduce random mutations ('ins','del','sub')

use Moose;

use CracTools::Const;
use CracTools::SimCT::Const;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;

# TODO create a parent class for all mutations generators
# that uses rates...

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

# Remove undefined settings
around 'BUILDARGS' => sub{
  my $orig  = shift;
  my $class = shift;
  my %args  = @_;
  for my $key ( keys %args ){
    delete $args{$key} unless defined $args{$key};
  }
  return $class->$orig(%args);
};

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
      $nb_sub-- if $self->genome_simulator->addMutation(
        CracTools::SimCT::Mutation::Substitution->new(
          chr => $chr,
          pos => int(rand($chr_length)),
          new_nuc => $CracTools::Const::NUCLEOTIDES->[int(rand(4))],
        ),
      );
    }

    # Generate Insertions
    while($nb_ins > 0) {
      $nb_ins-- if $self->genome_simulator->addMutation(
        CracTools::SimCT::Mutation::Insertion->new(
          chr => $chr,
          pos => int(rand($chr_length)),
          inserted_sequence => join("", 
            map{$CracTools::Const::NUCLEOTIDES->[int(rand(4))]} (1..int(rand($self->max_ins-1)+1))
          ), # This generates a random insertion
        ),
      );
    }

    # Generate Deletions
    while($nb_del > 0) {
      $nb_del-- if $self->genome_simulator->addMutation(
        CracTools::SimCT::Mutation::Deletion->new(
          chr => $chr,
          pos => int(rand($chr_length - $self->max_del)),
          length => int(rand($self->max_del-1)) + 1,
        ),
      );
    }
  }
}

1;

package CracTools::SimCT::MutationGenerator;  
# ABSTRACT: The mother class of mutation generators

use strict;
use warnings;

use Carp;

sub new {
  my $class = shift;
  my %args = @_;

  croak "Missing 'genome_simulator' argument" unless defined $args{genome_simulator};

  my $self = bless {
    genome_simulator => $args{genome_simulator},
  }, $class;

  $self->_init(\%args);

  return $self;
}

sub _init {}

sub genomeSimulator {
  my $self = shift;
  return $self->{genome_simulator};
}

sub generateMutations {
  croak "This class needs to be overrideen";
}

1;

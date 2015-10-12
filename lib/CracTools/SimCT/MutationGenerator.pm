package CracTools::SimCT::MutationGenerator;  
# ABSTRACT: The mother class of mutation generators

use Moose::Role;

use Carp;

has 'genome_simulator' => (
  is => 'rw',
  isa => 'CracTools::SimCT::GenomeSimulator',
  required => 1,
);

requires 'generateMutations';

1;

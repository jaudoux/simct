package CracTools::SimCT::MutationGenerator;  
# ABSTRACT: The mother class of mutation generators

use Moose::Role;

use Carp;

has 'genome_simulator';

requires 'generateMutations';

1;

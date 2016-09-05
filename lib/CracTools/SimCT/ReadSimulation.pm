package CracTools::SimCT::ReadSimulation;
# ABSTRACT: Base class for read simulations

use Moose::Role;

requires 'getGenomicIntervalsIterator', 'getSequenceIterator', 'getErrorsPosIterator', 'isPairedEnd';

has 'simulated_genome' => (
  is => 'ro',
  isa => 'CracTools::SimCT::SimulatedGenome',
  required => 1,
);

has 'genome_dir' => (
  is  => 'ro',
  isa => 'Str',
  required => 1,
);

has 'annotation_file' => (
  is  => 'ro',
  isa => 'Str',
);

has 'simulation_dir' => (
  is  => 'rw',
  isa => 'Str',
);

no Moose;
1;

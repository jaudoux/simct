package CracTools::SimCT::ReadSimulation;
# ABSTRACT: Base class for read simulations

use Moose::Role;

requires 'getGenomicIntervalsIterator', 'getSequenceIterator', 'getErrorsPosIterator', 'isPairedEnd';

no Moose;
1;

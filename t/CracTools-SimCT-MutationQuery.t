use strict;
use warnings;

use Test::More tests => 3;
use CracTools::SimCT::MutationQuery;
use CracTools::SimCT::Mutation::Substitution;

my $mut = CracTools::SimCT::Mutation::Substitution->new(
  chr => '1',
  start => 12,
  new_nuc => 'G'
);

my $mut_query = CracTools::SimCT::MutationQuery->new();
$mut_query->addMutation($mut);

{
  my $mutations = $mut_query->fetchByRegion('1',10,14);
  is(@{$mutations}, 1);
  is($mutations->[0], $mut);
}

{
  my $mutations = $mut_query->fetchByRegion('1',8,11);
  is(@{$mutations}, 0);
}

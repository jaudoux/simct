use strict;
use warnings;

use Test::More tests => 15;
use CracTools::SimCT::Mutation::Insertion;
use CracTools::SimCT::Mutation::Deletion;
use CracTools::SimCT::Mutation::Substitution;

my $chr_seq = "AGGCTGATCAGTCGATCAGTCAGTCAGCATCGAT";

# Test insertion
{
  my $ins = CracTools::SimCT::Mutation::Insertion->new(
    chr => "1",
    start => 10,
    inserted_sequence => "TATA",
  );

  is($ins->chr,"1");
  is($ins->start,10);
  is($ins->mutation_sequence,"TATA");
  is($ins->mutationLength,4);
  is($ins->referenceLength,0);
  is($ins->reference_sequence,"");
}

# Test deletion
{
  my $del = CracTools::SimCT::Mutation::Deletion->new(
    chr => "1",
    start => 10,
    length => 2,
  );
  is($del->mutation_sequence,"");
  is($del->mutationLength,0);
  is($del->referenceLength,2);
  is($del->reference_sequence,"NN");
  $del->setReferenceSequence("GT");
  is($del->reference_sequence,"GT");
}

# Test substitution
{
  my $sub = CracTools::SimCT::Mutation::Substitution->new(
    chr => "1",
    start => 10,
    new_nuc => 'C',
  );
  is($sub->mutation_sequence,'C');
  is($sub->mutationLength,1);
  is($sub->referenceLength,1);
  is($sub->reference_sequence,'N');
}

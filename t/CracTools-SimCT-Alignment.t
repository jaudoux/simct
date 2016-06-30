use strict;
use warnings;

use Test::More tests => 11;
use CracTools::SimCT::Alignment;
use CracTools::SimCT::Alignment::CigarElement;

sub newAlignment {
  my ($chr, $start, $query_length) = @_;
  return CracTools::SimCT::Alignment->new(chr => $chr, start => $start, query_length => $query_length);
}

sub newCigel {
  my ($op,$nb) = @_;
  return CracTools::SimCT::Alignment::CigarElement->new(op => $op, nb => $nb);
}

{
  my $alignment = newAlignment("1",10,40);
  is($alignment->chr, "1");
  is($alignment->start, 10);
  is($alignment->cigar, "40S");

  $alignment->appendCigarElement(newCigel("M",10));
  is($alignment->cigar,"10M30S");
  is($alignment->end, 19);

  $alignment->appendCigarElement(newCigel("N",20));
  $alignment->appendCigarElement(newCigel("M",10));
  $alignment->query_length(40);
  is($alignment->cigar,"10M20N10M20S");
  is($alignment->end, 49);
  is($alignment->query_length, 40);
  is($alignment->length, 40);

  $alignment->query_mapping_start(20);
  is($alignment->cigar,"20S10M20N10M");

  $alignment->query_mapping_start(10);
  is($alignment->cigar,"10S10M20N10M10S");
}

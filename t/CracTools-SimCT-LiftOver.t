use strict;
use warnings;

use Test::More tests => 18;
use CracTools::SimCT::LiftOver;

my $liftover = CracTools::SimCT::LiftOver->new();

#      [   10   ]  5 ] 4 ]    10   ]                       ]
# sim: AGTAGATGCATGTCAGTCA----------GTCAGCTAGCATGCATCGGGCTTA
# ref: AGTAGATGCA-----GTCACCCCCCCCCCGTCAGCTAGCATGCATCGGGCTTA

$liftover->addInterval("chr1",0,9,0);
$liftover->addInterval("chr1",15,18,-5); # 5bp insertion
$liftover->addInterval("chr1",19,25,+5); # 10bp deletion
$liftover->addInterval("chr1",26,40,0,"chr2"); # 10bp deletion

# Overlap of the insertion
{
  my $annot = $liftover->shiftAnnotation("chr1",2,16);
  is($annot->{chr},"chr1");
  is($annot->{start},2);
  is($annot->{end},11);
  is($annot->{strand},'+');
}

# Overlap of the deletion
{
  my $annot = $liftover->shiftAnnotation("chr1",16,20);
  is($annot->{chr},"chr1");
  is($annot->{start},11);
  is($annot->{end},25);
  is($annot->{strand},'+');
}

# Overlap of both insertion and deletion
{
  my $annot = $liftover->shiftAnnotation("chr1",2,20);
  is($annot->{chr},"chr1");
  is($annot->{start},2);
  is($annot->{end},25);
  is($annot->{strand},'+');
}

# Check alignement query
{
  my ($alignment,$chimeric_alignment) = $liftover->getAlignments("chr1",2,30);
  is($alignment->{chr},"chr1");
  is($alignment->{start},2);
  is($alignment->{cigar},"8M5I4M10D7M");
  is($chimeric_alignment->{chr},"chr2");
  is($chimeric_alignment->{start},26);
  is($chimeric_alignment->{cigar},"5M");
}

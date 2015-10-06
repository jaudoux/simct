use strict;
use warnings;

use Test::More tests => 12;
use CracTools::SimCT::LiftOver;

my $liftover = CracTools::SimCT::LiftOver->new();

#      [   10   ]  5 ] 4 ]    10   ]                       ]
# ref: AGTAGATGCATGTCAGTCA----------GTCAGCTAGCATGCATCGGGCTTA
# alt: AGTAGATGCA-----GTCACCCCCCCCCCGTCAGCTAGCATGCATCGGGCTTA

$liftover->addInterval("chr1",0,9,0);
$liftover->addInterval("chr1",15,18,-5); # 5bp deletion
$liftover->addInterval("chr1",19,25,+5); # 10bp insertion

# Overlap of the deletion
{
  my $annot = $liftover->shiftAnnotation("chr1",2,16);
  is($annot->{chr},"chr1");
  is($annot->{start},2);
  is($annot->{end},11);
  is($annot->{strand},'+');
}

# Overlap of the insertion
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

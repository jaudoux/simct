use strict;
use warnings;

use Test::More tests => 62;
#use Test::More tests => 4;
use CracTools::SimCT::LiftOver;
use CracTools::SimCT::GenomicInterval;

my $liftover = CracTools::SimCT::LiftOver->new();

sub newInterval {
  my ($chr,$start,$end,$strand) = @_;
  return CracTools::SimCT::GenomicInterval->new(
    chr   => $chr,
    start => $start,
    end   => $end,
    strand => defined $strand? $strand : '+',
  );
}

#      [   10   ]  5 ] 4 ]    10   ]                       ]
#                1111111111222222222233333333334444444444555
#      01234567890123456789012345678901234567890123456789012
#
#                111111111          1222222222233333333334444444444555
#      0123456789012345678          9012345678901234567890123456789012
# sim: AGTAGATGCATGTCAGTCA----------GTCAGCTAGCATGCATCGGGCTTA
# ref: AGTAGATGCA-----GTCACCCCCCCCCCGTCAGCTAGCATGCATCGGGCTTA
#      0123456789     0123456789012345678901234567890123456789012
#                     1111111111222222222233333333334444444444555

$liftover->addInterval("chr1",0,9,0);
$liftover->addInterval("chr1",15,18,-5); # 5bp insertion
$liftover->addInterval("chr1",19,25,+5); # 10bp deletion
$liftover->addInterval("chr1",26,30,-26,"chr2"); # chimeric alignment
$liftover->addInterval("chr1",31,40,40,"chr3",1); # chimeric alignment reversed (10 + 26)

# Overlap of the insertion
{
  my $annot = $liftover->shiftAnnotation(newInterval("chr1",2,16));
  is($annot->chr,"chr1");
  is($annot->start,2);
  is($annot->end,11);
  is($annot->strand,'+');
}

# Overlap of the deletion
{
  my $annot = $liftover->shiftAnnotation(newInterval("chr1",16,20));
  is($annot->chr,"chr1");
  is($annot->start,11);
  is($annot->end,25);
  is($annot->strand,'+');
}

# Included in the insertion
{
  my $annot = $liftover->shiftAnnotation(newInterval("chr1",10,14));
  is($annot, undef);
}

# Overlap of both insertion and deletion
{
  my $annot = $liftover->shiftAnnotation(newInterval("chr1",2,20));
  is($annot->chr,"chr1");
  is($annot->start,2);
  is($annot->end,25);
  is($annot->strand,'+');
}

{
  my $intervals = $liftover->shiftInterval(newInterval("chr1",2,30));
  is(scalar @{$intervals}, 4);

  is($intervals->[0]->chr,"chr1");
  is($intervals->[0]->start,2);
  is($intervals->[0]->end,9);
  is($intervals->[0]->reference_interval->chr,"chr1");
  is($intervals->[0]->reference_interval->start,2);
  is($intervals->[0]->reference_interval->end,9);

  is($intervals->[1]->chr,"chr1");
  is($intervals->[1]->start,10);
  is($intervals->[1]->end,13);
  is($intervals->[1]->reference_interval->chr,"chr1");
  is($intervals->[1]->reference_interval->start,15);
  is($intervals->[1]->reference_interval->end,18);

  is($intervals->[2]->chr,"chr1");
  is($intervals->[2]->start,24);
  is($intervals->[2]->end,30);
  is($intervals->[2]->reference_interval->chr,"chr1");
  is($intervals->[2]->reference_interval->start,19);
  is($intervals->[2]->reference_interval->end,25);

  is($intervals->[3]->chr,"chr2");
  is($intervals->[3]->start,0);
  is($intervals->[3]->end,4);
  is($intervals->[3]->reference_interval->chr,"chr1");
  is($intervals->[3]->reference_interval->start,26);
  is($intervals->[3]->reference_interval->end,30);
}

{
  my $intervals = $liftover->shiftInterval(newInterval("chr1",31,35,'-'));
  is(scalar @{$intervals}, 1);
  is($intervals->[0]->chr,"chr3");
  is($intervals->[0]->strand,'+');
  is($intervals->[0]->start,5);
  is($intervals->[0]->end,9);
}

# Check alignement query
{
  my ($alignment,$chimeric_alignment) = $liftover->getAlignments(newInterval("chr1",2,30));
  is($alignment->chr,"chr1");
  is($alignment->start,2);
  is($alignment->cigar,"8M5I4M10D7M5S");
  is($chimeric_alignment->chr,"chr2");
  is($chimeric_alignment->start,0);
  is($chimeric_alignment->cigar,"24S5M");
}

# Check reverse chimeric alignement query
{
  my ($chimeric_alignment,$reverse_chimeric_alignment) = $liftover->getAlignments(newInterval("chr1",28,36,'-'));
  is($chimeric_alignment->chr,"chr2");
  is($chimeric_alignment->start,2);
  is($chimeric_alignment->cigar,"3M6S");
  is($chimeric_alignment->strand,'-');
  is($reverse_chimeric_alignment->chr,"chr3");
  is($reverse_chimeric_alignment->strand, '+');
  is($reverse_chimeric_alignment->start,4);
}

# Check spliced alignement query
{
  my ($alignment,$chimeric_alignment) = $liftover->getSplicedAlignments(
    newInterval("chr1",2,10),
    newInterval("chr1",20,30),
  );
  is($alignment->chr,"chr1");
  is($alignment->start,2);
  is($alignment->cigar,"8M1I15N6M5S");
  is($chimeric_alignment->chr,"chr2");
  is($chimeric_alignment->start,0);
  is($chimeric_alignment->cigar,"15S5M");
}

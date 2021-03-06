use strict;
use warnings;

use Test::More tests => 145;
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
$liftover->addInterval("chr1",41,100,0); # chimeric alignment reversed (10 + 26)
$liftover->addInterval("chr2",0,100,0);



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

# Check spliced alignement on forward strand
{
  my @a = $liftover->getSplicedAlignments(newInterval("chr2",1,5,'+'),newInterval("chr2",10,20,'+'));
  is($a[0]->chr,"chr2");
  is($a[0]->start,1);
  is($a[0]->cigar,"5M4N11M");
  is($a[0]->strand,'+');
}

# Check false-spliced alignement overlapping an insertion
{
  my @a = $liftover->getSplicedAlignments(newInterval("chr1",1,11,'+'),newInterval("chr1",15,18,'+'));
  is(scalar @a, 1);
  is($a[0]->chr, "chr1");
  is($a[0]->start, 1);
  is($a[0]->strand,'+');
  is($a[0]->cigar, "9M2I4M");
}

# Check true-spliced alignement overlapping an insertion
{
  my @a = $liftover->getSplicedAlignments(newInterval("chr1",1,11,'+'),newInterval("chr1",19,22,'+'));
  is(scalar @a, 1);
  is($a[0]->chr, "chr1");
  is($a[0]->start, 1);
  is($a[0]->strand,'+');
  is($a[0]->cigar, "9M2I14N4M");
}


# Check multi spliced alignement on forward strand
{
  my @a = $liftover->getSplicedAlignments(newInterval("chr2",1,5,'+'),newInterval("chr2",8,12,'+'),newInterval("chr2",15,20,'+'));
  is($a[0]->chr,"chr2");
  is($a[0]->start,1);
  is($a[0]->cigar,"5M2N5M2N6M");
  is($a[0]->strand,'+');
}

# Check spliced alignement on reverse strand
{
  my @a = $liftover->getSplicedAlignments(newInterval("chr2",1,5,'-'),newInterval("chr2",10,20,'-'));
  is($a[0]->chr,"chr2");
  is($a[0]->start,1);
  is($a[0]->cigar,"5M4N11M");
  is($a[0]->strand,'-');
}

# Simple fusion simulation
$liftover->addInterval("fusion_1",0,49,0,"chr1"); # Start at chr1 pos 0
$liftover->addInterval("fusion_1",50,100,-50,"chr2"); # Start at chr2 pos 0

# Check chimeric alignment query
{
  my @a = $liftover->getAlignments(newInterval("fusion_1",25,74,'+'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '+');
  is($a[0]->cigar, '25M25S');
  is($a[1]->start, 0);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '+');
  is($a[1]->cigar, '25S25M');
}

# Check spliced chimeric alignment query
{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_1",25,74,'+'), newInterval("fusion_1",90,100,'+'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '+');
  is($a[0]->cigar, '25M36S');
  is($a[1]->start, 0);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '+');
  is($a[1]->cigar, '25S25M15N11M');
}

# Check reverse chimeric alignment query
{
  my @a = $liftover->getAlignments(newInterval("fusion_1",25,74,'-'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '-');
  is($a[0]->cigar, '25M25S');
  is($a[1]->start, 0);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '-');
  is($a[1]->cigar, '25S25M');
}

# Check reverse spliced chimeric alignment query
{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_1",25,74,'-'), newInterval("fusion_1",90,100,'-'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '-');
  is($a[0]->cigar, '25M36S');
  is($a[1]->start, 0);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '-');
  is($a[1]->cigar, '25S25M15N11M');
}

# Forward-Reverse fusion simulation
$liftover->addInterval("fusion_2",0,49,0,"chr1"); # Start at chr1 pos 0
$liftover->addInterval("fusion_2",50,100,100,"chr2",1); # Start at chr2 pos 0

# Check chimeric alignment query
{
  my @a = $liftover->getAlignments(newInterval("fusion_2",25,74,'+'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '+');
  is($a[0]->cigar, '25M25S');
  is($a[1]->start, 26);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '-');
  is($a[1]->cigar, '25S25M');
}

# Check spliced chimeric alignment query
{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_2",25,74,'+'), newInterval("fusion_2",90,100,'+'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '+');
  is($a[0]->cigar, '25M36S');
  is($a[1]->start, 0);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '-');
  is($a[1]->cigar, '11M15N25M25S');
}

# Check reverse chimeric alignment query
{
  my @a = $liftover->getAlignments(newInterval("fusion_2",25,74,'-'));
  is($a[0]->start, 25);
  is($a[0]->chr, "chr1");
  is($a[0]->strand, '-');
  is($a[0]->cigar, '25M25S');
  is($a[1]->start, 26);
  is($a[1]->chr, "chr2");
  is($a[1]->strand, '+');
  is($a[1]->cigar, '25S25M');
}

# class3 fusion simulation
$liftover->addInterval("fusion_3",0,49,49,"chr1",1); # Start at chr1 pos 50
$liftover->addInterval("fusion_3",50,100,250,"chr1",1); # Start at chr1 pos 0
{
  my @a = $liftover->getAlignments(newInterval("fusion_3",25,74,'+'));
  is($a[0]->start, 0);
  is($a[0]->strand, '-');
  is($a[0]->cigar, '25M25S');
  is($a[1]->start, 176);
  is($a[1]->strand, '-');
  is($a[1]->cigar, '25S25M');
}

{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_3",25,40,'+'), newInterval("fusion_3",60,75,'+'));
  is(scalar @a, 2);
}

# reverse splice fusion simulation
$liftover->addInterval("fusion_4",0,49,250,"chr1",1); # Start at chr1 pos 50
$liftover->addInterval("fusion_4",50,100,100,"chr1",1); # Start at chr1 pos 0
{
  my @a = $liftover->getAlignments(newInterval("fusion_4",25,74,'-'));
  is(scalar @a, 1);
  is($a[0]->cigar, '25M150N25M');
  is($a[0]->start, 26); # 100 - 74
}

{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_4",25,40,'+'), newInterval("fusion_4",65,75,'+'));
  is(scalar @a, 1);
  # use Data::Dumper;
  # print STDERR Dumper(\@a);
  is($a[0]->start, 25);
  is($a[0]->cigar, "11M174N16M");
}

{
  my @a = $liftover->getSplicedAlignments(newInterval("fusion_4",0,10,'-'), newInterval("fusion_4",25,70,'-'));
  is(scalar @a, 1);
  is($a[0]->start, 30);
  is($a[0]->cigar, "21M150N25M14N11M");
}



# Class 3 fusion
$liftover->addInterval("fusion_5",0,49,250,"chr1",1); # Start at chr1 pos 0
$liftover->addInterval("fusion_5",50,100,-50,"chr1"); # Start at chr2 pos 0
{
  my @a = $liftover->getAlignments(newInterval("fusion_5",25,74,'-'));
  is(scalar @a, 2);
}

# Class 3 fusion
$liftover->addInterval("fusion_6",0,49,49,"chr1",1); # Start at chr1 pos 0
$liftover->addInterval("fusion_6",50,100,250,"chr1"); # Start at chr2 pos 0
{
  my @a = $liftover->getAlignments(newInterval("fusion_6",25,74,'-'));
  is(scalar @a, 2);
}

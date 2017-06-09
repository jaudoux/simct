#! /usr/bin/perl
#
use strict;
use warnings;

my $fh;

while(<>) {
  if($_ =~ /^>(\S+)/) {
    open($fh, '>', "chr".$1.".fa") or die "Cannot open file";
  }
  print $fh $_;
}

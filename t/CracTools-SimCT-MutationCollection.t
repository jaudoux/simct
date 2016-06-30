use strict;
use warnings;

use Test::More tests => 14;
use CracTools::SimCT::MutationCollection;
use Inline::Files 0.68;
use File::Temp;

{
  # Create a temp GTF annotation file
  my $vcf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<VCF>) {print $vcf_file $_;}
  close $vcf_file;

  my $mutation_collection = CracTools::SimCT::MutationCollection->new();

  $mutation_collection->loadVCF($vcf_file);

  is(scalar @{$mutation_collection->substitutions}, 4);
  is(scalar @{$mutation_collection->insertions}, 1);
  is(scalar @{$mutation_collection->deletions}, 1);

  my $first_snv = $mutation_collection->substitutions->[0];
  is($first_snv->chr, '20');
  is($first_snv->start, 14370);
  # For now, the reference sequence is setted to 'N' by default
  #is($first_snv->reference_sequence, 'G');
  is($first_snv->mutation_sequence, 'A');

  my $first_ins = $mutation_collection->insertions->[0];
  is($first_ins->chr, '20');
  is($first_ins->start, 1234571);
  is($first_ins->mutation_sequence, 'T');
  is($first_ins->referenceLength, 0);

  my $first_del = $mutation_collection->deletions->[0];
  is($first_del->chr, '20');
  is($first_del->start, 1234568);
  is($first_del->mutation_sequence, '');
  is($first_del->referenceLength, 3);
}

__VCF__
##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA2        NA00002        NA00003
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4
20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	1234567	microsat1	GTCT	G,GTCTT	50	PASS	NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3

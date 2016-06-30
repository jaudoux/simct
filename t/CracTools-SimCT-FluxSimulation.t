package TestSimulation;
use Moose;
use CracTools::SimCT::GenomicInterval;
with 'CracTools::SimCT::ReadSimulation';

sub getGenomicIntervalsIterator {
  my @intervals = (
    CracTools::SimCT::GenomicInterval->new(chr => 1, start => 1, end => 12),
    CracTools::SimCT::GenomicInterval->new(chr => 1, start => 2, end => 13),
  );
  return sub {
    my $interval = shift @intervals;
    return $interval;
  }
}

sub getSequenceIterator {
  my @sequences = (
    { name => 'read1', seq => 'GATGCTAGTCGA', qual => '############'},
    { name => 'read2', seq => 'ATGCTAGTCGAT', qual => '############'},
  );
  return sub {
    my $seq = shift @sequences;
    #use Data::Dumper;
    #print STDERR Dumper($seq);
    return $seq;
  }
}

sub getErrorsPosIterator {
  return sub {
    return ();
  }
}

sub isPairedEnd {
  return 0;
}

no Moose;

#------------------------------------------------------------------------------#
package TestSimulator;
use Moose;
use TestSimulation;
with 'CracTools::SimCT::ReadSimulator';

has '+simulator_name' => (
  default => 'test_simulator',
);

sub _generateSimulation {
  return TestSimulation->new(@_);
}

no Moose;

#------------------------------------------------------------------------------#
package test;

use strict;
use warnings;

use Test::More tests => 5;

use CracTools::Utils;
use CracTools::SimCT::Annotations;
use CracTools::SimCT::Genome;
use CracTools::SimCT::GenomeSimulator;
use CracTools::SimCT::Mutation::Substitution;
use CracTools::SimCT::ReadSimulator::FluxSimulator;

use File::Spec;

use Inline::Files 0.68;
use File::Temp;

{
  # Create a temp GTF annotation file
  my $gtf_file = new File::Temp( SUFFIX => '.gtf', UNLINK => 1);
  while(<ANNOTATIONS>) {print $gtf_file $_;}
  close $gtf_file;

  # Create a temp genome directory and print the FASTA files
  my $chr1_file = new File::Temp( SUFFIX => '.fa', UNLINK => 1);
  while(<CHR1>) {print $chr1_file $_;}
  close $chr1_file;

  # Create the genome and the genome simulator
  my $genome = CracTools::SimCT::Genome->new(
    reference_sequence_files => { "1" => $chr1_file->filename }
  );
  my $genome_sim = CracTools::SimCT::GenomeSimulator->new(genome => $genome);

  # Load annotations
  my $annotations = CracTools::SimCT::Annotations->new();
  $annotations->loadGTF($gtf_file);

  # Add one mutation
  $genome_sim->addMutation(
    CracTools::SimCT::Mutation::Substitution->new(
      chr => 1, start => 13, new_nuc => 'G'),
  );

  # Simulate the genome
  my $simulated_genome_dir = File::Temp->newdir();
  my ($simulated_annotations, $simulated_genome) = $genome_sim->generateGenome(
    genome_dir => $simulated_genome_dir,
    annotations => $annotations,
  );

  # Run the simulation with the "test" simulator
  my $test_simulator = TestSimulator->new();
  my $simulation_dir = File::Temp->newdir();
  my $test_simulation = $test_simulator->runSimulation(
    simulated_genome => $simulated_genome,
    genome_dir => $simulated_genome_dir,
    simulation_dir => $simulation_dir,
  );

  is($test_simulator->simulator_name,'test_simulator');
  my $simulation_FASTQ_file = File::Spec->catfile($simulation_dir,"reads.fastq.gz");
  my $fastq_it = CracTools::Utils::seqFileIterator($simulation_FASTQ_file,'fastq');
  my $read1 = $fastq_it->();
  is($read1->{name},"0:1,1,+,12M");
  is($read1->{seq},"GATGCTAGTCGA");

  my $read2 = $fastq_it->();
  is($read2->{seq},"ATGCTAGTCGAT");
  is($read2->{name},"1:1,2,+,12M");

  # Run the simulation with flux
  # my $flux_simulator = CracTools::SimCT::ReadSimulator::FluxSimulator->new();
  # my $simulation_dir = File::Temp->newdir();
  # my $flux_simulation = $flux_simulator->runSimulation(
  #   simulated_genome => $simulated_genome,
  #   genome_dir => $simulated_genome_dir,
  #   simulation_dir => $simulation_dir,
  #   annotation_file => $simulated_annotations,
  #   flux_parameters => {
  #     READ_NUMBER       => 10,
  #     NB_MOLECULES      => 10,
  #     READ_LENGTH       => 10,
  #   },
  # );
}


__CHR1__
>1
GATGCTAGTCGATGATATTTTATAGCTAGCATGGGATCATCGATG
__ANNOTATIONS__
1	StringTie	transcript	10	25	1000	+	.	gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "105.907524"; FPKM "128.387177";
1	StringTie	exon	10	15	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "1"; cov "42.406654";
1	StringTie	exon	20	25	1000	+	.	gene_id "geneA"; transcript_id "geneA.1"; exon_number "2"; cov "161.272552";

package CracTools::SimCT::ReadSimulator::FluxSimulator;
# ABSTRACT: Base class for read simulators

use Moose;

use Carp;
use File::Spec;
use CracTools::Utils;
use CracTools::SimCT::Const;
use CracTools::SimCT::ReadSimulation::FluxSimulator;

with 'CracTools::SimCT::ReadSimulator';

has '+simulator_name' => (
  default => 'FluxSimulator',
);

has 'flux_binary' => (
  is => 'rw',
  isa => 'Str',
  default => $CracTools::SimCT::Const::FLUX_BINARY,
);

around '_postProcessSimulation' => sub {
  my $orig = shift;
  my $self = shift;
  my $simulation = shift;
  my $simulation_dir = $simulation->simulation_dir;

  # Print gene_counts and transcript_counts files
  my $transcript_counts_output = File::Spec->catfile($simulation_dir,"transcript-counts.tsv.gz");
  my $transcript_counts_fh     = CracTools::Utils::getWritingFileHandle($transcript_counts_output);
  my $gene_counts_output       = File::Spec->catfile($simulation_dir,"gene-counts.tsv.gz");
  my $gene_counts_fh           = CracTools::Utils::getWritingFileHandle($gene_counts_output);
  my $flux_profile_fh          = CracTools::Utils::getReadingFileHandle($simulation->{profile_file});
  my $gtf_it                   = CracTools::Utils::gffFileIterator($simulation->annotation_file,'gtf');

  # 1. Create a transcript_id => gene conversion hash
  my %conversion_hash;
  while(my $annot = $gtf_it->()) {
    next if $annot->{feature} ne 'exon';
    my $gene_id = $annot->{attributes}->{gene_id};
    my $transcript_id = $annot->{attributes}->{transcript_id};
    next if !defined $gene_id || !defined $transcript_id;
    $conversion_hash{$transcript_id} = $gene_id;
  }

  # 2. Merge transrcipt counts into gene counts and print transcript counts
  print $transcript_counts_fh join("\t","feature","truth"),"\n";
  tie my %gene_counts, 'Tie::RefHash';
  while(<$flux_profile_fh>) {
    my @fields = split("\t",$_);
    # TODO HANDLE FUSION TRANSCRIPTS!!!!
    next if !defined $conversion_hash{$fields[1]};
    $gene_counts{$conversion_hash{$fields[1]}} += $fields[9];
    print $transcript_counts_fh join("\t",$fields[1],$fields[9]),"\n";
  }

  # 3. Print gene counts to the corresponding output file
  print $gene_counts_fh join("\t","feature","truth"),"\n";
  foreach my $gene (keys %gene_counts) {
    print $gene_counts_fh join("\t",$gene,$gene_counts{$gene}),"\n";
  }

  return $self->$orig($simulation);
};

sub _generateSimulation {
  my $self            = shift;
  my %args            = @_;

  my $simulated_genome = $args{simulated_genome};
  my $output_dir       = $args{output_dir};
  my $genome_dir       = $args{genome_dir};
  my $annotation_file  = $args{annotation_file};
  my $flux_parameters  = $args{flux_parameters};

  $flux_parameters = {} unless defined $flux_parameters;

  croak "Mission 'genome_dir' argument" unless defined $genome_dir;
  croak "Mission 'output_dir' argument" unless defined $output_dir;
  croak "Mission 'annotation_file' argument" unless defined $annotation_file;

  $flux_parameters->{FASTA}       = 'YES'; # Force sequence output
  $flux_parameters->{ERR_FILE}    = 76 unless defined $flux_parameters->{ERR_FILE}; # Force sequence output
  $flux_parameters->{UNIQUE_IDS}  = 'YES' if defined $flux_parameters->{PAIRED_END} && $flux_parameters->{PAIRED_END} eq 'YES';

  my $flux_output_name  = $CracTools::SimCT::Const::FLUX_OUTPUT_BASENAME;
  my $parameter_file    = File::Spec->catfile($output_dir,"$flux_output_name.par");
  my $fastq_file        = File::Spec->catfile($output_dir,"$flux_output_name.fastq");
  my $library_file      = File::Spec->catfile($output_dir,"$flux_output_name.lib");
  my $profile_file      = File::Spec->catfile($output_dir,"$flux_output_name.pro");
  my $sequencing_file   = File::Spec->catfile($output_dir,"$flux_output_name.bed");

  # First we need to create a parameter file, as expected for Flux
  my $parameter_fh = CracTools::Utils::getWritingFileHandle($parameter_file);
  print $parameter_fh "GEN_DIR\t$genome_dir\n";
  print $parameter_fh "REF_FILE_NAME\t$annotation_file\n";
  foreach my $parameter (keys %{$flux_parameters}) {
    print $parameter_fh "$parameter\t".$flux_parameters->{$parameter}."\n";
  }
  close($parameter_fh);

  # Remove sorted annotation that may exists
  unlink File::Spec->catfile($output_dir,"annotations_sorted.gtf");

  # Then run flux
  my $command = $self->flux_binary." -p $parameter_file --force -x -l -s";
  system($command);

  return CracTools::SimCT::ReadSimulation::FluxSimulator->new(
    annotation_file  => $annotation_file,
    simulated_genome => $simulated_genome,
    genome_dir       => $genome_dir,
    flux_parameters  => $flux_parameters,
    profile_file     => $profile_file,
    parameter_file   => $parameter_file,
    library_file     => $library_file,
    sequencing_file  => $sequencing_file,
    fastq_file       => $fastq_file,
  );
}

no Moose;
__PACKAGE__->meta->make_immutable;

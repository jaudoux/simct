package CracTools::SimCT::FluxWrapper;
# ABSTRACT: Wrapper around FluxSimulator binary tool

use strict;
use warnings;

use Carp;
use File::Spec;
use CracTools::Utils;
use CracTools::SimCT::Const;

sub new {
  my $class = shift;
  my %args = @_;

  my $flux_binary = $args{flux_binary};
  $flux_binary = $CracTools::SimCT::Const::FLUX_BINARY unless defined $flux_binary;

  my $self = bless {
    flux_binary => $flux_binary,
  }, $class;

  #$self->_init();

  return $self;

}

#sub _init {
#
#}

sub fluxBinary {
  my $self = shift;
  return $self->{flux_binary};
}

# Given a genome and an annotation file object, generate a set of read files and other available files :
sub generateSimulation {
  my $self            = shift;
  my %args            = @_;

  my $genome_dir      = $args{genome_dir};
  my $annotation_file = $args{annotation_file};
  my $output_dir      = $args{output_dir};
  my $flux_parameters = $args{flux_parameters};

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
  my $command = $self->fluxBinary." -p $parameter_file --force -x -l -s";
  system($command);

  return {
    flux_parameters => $flux_parameters,
    profile_file    => $profile_file,
    parameter_file  => $parameter_file,
    library_file    => $library_file,
    sequencing_file => $sequencing_file,
    fastq_file      => $fastq_file,
    #error_file      => $error_file,
  }
}

sub getErrorsPos($) {
  my $seq = shift;
  my @pos;
  while ($seq =~ /[atgc]/g) {
    push @pos,$-[0];
  }
  return \@pos;
}

1;

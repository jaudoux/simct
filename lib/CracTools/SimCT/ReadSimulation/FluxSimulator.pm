package CracTools::SimCT::ReadSimulation::FluxSimulator;
# ABSTRACT: FluxSimulator simulation

use Moose;

with 'CracTools::SimCT::ReadSimulation';

has '+annotation_file' => (
  required => 1,
);

has 'flux_parameters' => (
  is => 'ro',
  isa => 'HashRef',
);

has 'profile_file' => (
  is => 'ro',
  isa => 'Str',
);

has 'parameter_file' => (
  is => 'ro',
  isa => 'Str',
);

has 'library_file' => (
  is => 'ro',
  isa => 'Str',
);

has 'sequencing_file' => (
  is => 'ro',
  isa => 'Str',
);

has 'fastq_file' => (
  is => 'ro',
  isa => 'Str',
);


sub getGenomicIntervalsIterator {
  my $self = shift;
  my $bed_it = CracTools::Utils::bedFileIterator($self->sequencing_file);
  return sub {
    my $bed_line = $bed_it->();
    if(defined $bed_line) {
      my @intervals;
      foreach my $block (@{$bed_line->{blocks}}) {
        # TODO we could do the "real strand" conversion here, nop?
        push @intervals, CracTools::SimCT::GenomicInterval->new(
          chr => $bed_line->{chr},
          start => $block->{ref_start},
          end => $block->{ref_end} - 1, # because BED intervals are half-open
          strand => $bed_line->{strand},
        );
      }
      return @intervals;
    } else {
      return undef;
    }
  }
}

sub getSequenceIterator {
  my $self = shift;
  my $seq_it = CracTools::Utils::seqFileIterator($self->fastq_file,'fastq');
  my $read_id = 0;
  return sub {
    my $read = $seq_it->();
    if(defined $read) {
      if($self->isPairedEnd && $read_id % 2 == 0) {
        $read->{reversed} = 1;
      } else {
        $read->{reversed} = undef;
      }
      $read->{read_id} = $read_id;
      $read_id++;
    }
    return $read;
  }
}

sub getErrorsPosIterator {
  my $self = shift;
  my $seq_it = CracTools::Utils::seqFileIterator($self->fastq_file,'fastq');
  return sub {
    my $read = $seq_it->();
    return defined $read? _getErrorsPos($read->{seq}) : undef;
  }
}

sub isPairedEnd {
  my $self = shift;
  return $self->flux_parameters->{PAIRED_END} eq 'YES'? 1 : 0;
}

sub _getErrorsPos($) {
  my $seq = shift;
  my @pos;
  while ($seq =~ /[atgc]/g) {
    push @pos,$-[0];
  }
  return @pos;
}

no Moose;
__PACKAGE__->meta->make_immutable;

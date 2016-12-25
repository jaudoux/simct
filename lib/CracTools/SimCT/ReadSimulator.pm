package CracTools::SimCT::ReadSimulator;
# ABSTRACT: Base class for read simulators

use Moose::Role;

use File::Path;
use File::Spec;
use Tie::RefHash;

use CracTools::Utils;
use CracTools::Output;
use CracTools::SimCT::Const;

has 'simulator_name' => (
  is => 'ro',
  isa => 'Str',
  init_arg => undef,
);

has 'uniq_ids' => (
  is => 'rw',
  isa => 'Bool',
  default => $CracTools::SimCT::Const::UNIQ_IDS,
);

has 'disable_error_encoding' => (
  is => 'rw',
  isa => 'Bool',
  default => $CracTools::SimCT::Const::DISABLE_ERROR_ENCODING,
);

# Create a new simulation object, this is overloaded by the base classes
# to create the right simulation object
requires '_generateSimulation';

# TODO we should be able to handle the simulation of multiple simulated genome
# that derives from the same initial genome.
sub runSimulation {
  my $self = shift;
  my %args = @_;
  my $simulated_genome = $args{simulated_genome};
  my $genome_dir       = $args{genome_dir};
  my $simulation_dir   = $args{simulation_dir};

  # Set the directory for the simulator
  $args{output_dir} = File::Spec->catfile($simulation_dir,$self->simulator_name);
  File::Path::make_path($args{output_dir});

  # Create the simulation
  my $simulation = $self->_generateSimulation(%args);
  $simulation->simulation_dir($simulation_dir);

  # A hook to handle the post-processing of the simulation
  return $self->_postProcessSimulation($simulation);
}

sub _postProcessSimulation {
  my $self = shift;
  my $simulation = shift;

  my $simulation_dir   = $simulation->simulation_dir;
  my $simulated_genome = $simulation->simulated_genome;

  # Open output files for reads
  my $fastq1_output = $simulation->isPairedEnd?
    File::Spec->catfile($simulation_dir,"reads_1.fastq.gz") :
    File::Spec->catfile($simulation_dir,"reads.fastq.gz");
  my $fastq2_output = $simulation->isPairedEnd?
    File::Spec->catfile($simulation_dir,"reads_2.fastq.gz") : undef;

  my $fastq1_output_fh = CracTools::Utils::getWritingFileHandle($fastq1_output);
  my $fastq2_output_fh = $simulation->isPairedEnd?
    CracTools::Utils::getWritingFileHandle($fastq2_output) : undef;

  # Get iterators for sequences, genomic intervals and error positions
  my $seq_it       = $simulation->getSequenceIterator();
  my $intervals_it = $simulation->getGenomicIntervalsIterator();
  my $errors_it    = $simulation->getErrorsPosIterator();

  # Set counters and loop vars
  my $nb_reads  = 0;
  my $nb_errors = 0;
  my @paired_alignments;
  my @paired_errors_pos;
  my $paired_seq;
  my $paired_qual;

  # Create a hash to store mutations, with their associated reads
  my %splices;
  my %chimeras;
  tie my %mutations, 'Tie::RefHash';

  # Loop over reads to write them in the output file with their alignements
  # and error positions as the read_name
  while (my $read = $seq_it->()) {

    # Set the read id
    my $read_id = $nb_reads;
    if($simulation->isPairedEnd) {
      $read_id = ($nb_reads % 2 == 0)? $nb_reads / 2 : ($nb_reads - 1) / 2;
    }

    # TODO add some hooks for subclasses
    my @intervals   = $intervals_it->();
    my @alignments  = $simulated_genome->liftover->getSplicedAlignments(@intervals);
    my @errors_pos  = $errors_it->();

    my @overlapping_mutations = ();

    # Add splices and/or chimeras (based on alignments)
    my $prev_alignment;
    foreach my $alignment (@alignments) {
      my $start = $alignment->start;
      my $real_strand = $alignment->strand;
      if($read->{reversed}) {
        $real_strand = $alignment->strand eq '+'? '-' : '+';
      }
      foreach my $cigel ($alignment->allCigarElements) {
        if($cigel->op eq 'N') {
          # Add splice
          if($cigel->nb <= $CracTools::SimCT::Const::MAX_SPLICE_LENGTH) {
            my $splice_key = join("@", $alignment->chr, $start,
              $start + $cigel->nb, $real_strand);
            push @{$splices{$splice_key}},$read_id;
          # Splice is larger than MAX_SPLICE_LENGTH, then it is a class2 chimera
          } else {
            my $chim_key = join('@', $alignment->chr, $start, $real_strand,
              $alignment->chr, $start + $cigel->nb, $real_strand);
            push @{$chimeras{$chim_key}},$read_id;
          }
          $start += $cigel->nb;
        # If this cigar op is reference based, we update
        # the start position
        } elsif($cigel->op =~ /^[MDX=]$/) {
          # Add overlapping mutations
          push @overlapping_mutations, $simulated_genome->mutation_query->getOverlappingMutations(
            CracTools::SimCT::GenomicInterval->new(
              chr   => $alignment->chr,
              start => $start,
              end   => $start + $cigel->nb - 1
          ));
          $start += $cigel->nb;
        }
      }
      if(defined $prev_alignment) {
        my $chim_key;

        $chim_key = join("@",
          $prev_alignment->chr,
          #$prev_alignment->strand eq '-' ? $prev_alignment->start : $prev_alignment->end,
          $real_strand ne $prev_alignment->strand ? $prev_alignment->start : $prev_alignment->end,
          $real_strand ne $prev_alignment->strand ? $prev_alignment->strand : $real_strand,
          $alignment->chr,
          #$alignment->strand eq '-' ? $alignment->end : $alignment->start,
          $real_strand ne $alignment->strand ? $alignment->end : $alignment->start,
          $real_strand ne $alignment->strand ? $alignment->strand : $real_strand);

        # if($read->{reversed}) {
        #   $chim_key = join("@",
        #     $alignment->chr, $alignment->start, $alignment->strand eq '+'? '-' : '+',
        #     $prev_alignment->chr, $prev_alignment->end, $prev_alignment->strand eq '+'? '-' : '+');
        # } else {
        #   $chim_key = join("@",
        #     $prev_alignment->chr, $prev_alignment->end, $prev_alignment->strand,
        #     $alignment->chr, $alignment->start, $alignment->strand);
        # }
        push @{$chimeras{$chim_key}},$read_id;
      }
      $prev_alignment = $alignment;
    }

    # Add mutation for this read to the main hash
    map { push @{$mutations{$_}}, $read_id } @overlapping_mutations;

    $nb_errors += scalar @errors_pos;

    # Handle uniq_ids
    if($simulation->isPairedEnd && $self->uniq_ids) {
      # This is the first pair, we only record the alignments and errors
      if($nb_reads % 2 == 0) {
        @paired_alignments = @alignments;
        @paired_errors_pos = @errors_pos;
        $paired_seq        = $read->{seq};
        $paired_qual       = $read->{qual};

      # This is the second pair, we merge alignments and errors and print
      # both reads to their respective files with the same read_name
      } else {

        # Append alignments and errors from the pair to the current read
        unshift @alignments, @paired_alignments;
        unshift @errors_pos, @paired_errors_pos;

        # Special lossy treatment for errors
        my @sorted_uniq_errors = do { my %seen; grep { !$seen{$_}++ } sort @errors_pos };

        # Encode the read name
        my $read_name = $self->disable_error_encoding?
          _getReadName($read_id, \@alignments) :
          _getReadName($read_id, \@alignments, \@sorted_uniq_errors);

        # Write both pairs to their respective file
        CracTools::Utils::writeSeq($fastq1_output_fh, 'fastq', $read_name."/1", $paired_seq,$paired_qual);
        CracTools::Utils::writeSeq($fastq2_output_fh, 'fastq', $read_name."/2", $read->{seq},$read->{qual});
      }

    # Default output (no uniq_ids)
    } else {
      # Get the right filehandle
      my $fh = ($nb_reads % 2 == 0) || !$simulation->isPairedEnd? $fastq1_output_fh : $fastq2_output_fh;
      # Write the seq to the output file
      CracTools::Utils::writeSeq($fh,'fastq',
        $self->disable_error_encoding?
          _getReadName($read_id, \@alignments) :
          _getReadName($read_id, \@alignments, \@errors_pos),
        $read->{seq},
        $read->{qual}
      );
    }
    $nb_reads++;
  }

  # Print Mutations (if any)
  if(keys %mutations) {
    my $mutations_output = File::Spec->catfile($simulation_dir,"mutations.vcf.gz");
    my $mutations_fh     = CracTools::Utils::getWritingFileHandle($mutations_output);
    foreach my $mut (sort { $a->chr cmp $b->chr || $a->start <=> $b->start } keys %mutations) {
      # Get raw VCF record
      my $vcf_line = $mut->getVCFRecord;
      # Update information with read names and depth
      $vcf_line->{id} = join(',',@{$mutations{$mut}});
      $vcf_line->{info}->{DP} = [scalar @{$mutations{$mut}}];
      CracTools::SimCT::Utils::printVCFLine($mutations_fh,$vcf_line);
    }
  }

  # Print Chimeras
  # TODO sort chimeras by pos
  if(keys %chimeras) {
    my $chimeras_output = File::Spec->catfile($simulation_dir,"chimeras.tsv.gz");
    my $chimeras_fh     = CracTools::Utils::getWritingFileHandle($chimeras_output);
  # TODO Sort chimeras by start position
    foreach my $chimera (keys %chimeras) {
      my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = split("@",$chimera);
      print $chimeras_fh join("\t",
        $chr1,$pos1,$strand1,$chr2,$pos2,$strand2,
        join(':',@{$chimeras{$chimera}}),
        scalar @{$chimeras{$chimera}},
      ),"\n";
    }
  }

  # Print Splices
  # TODO sort splices by pos
  if(keys %splices) {
    my $splices_output = File::Spec->catfile($simulation_dir,"splices.bed.gz");
    my $splices_fh     = CracTools::Utils::getWritingFileHandle($splices_output);
  # TODO Sort splices by start position
    foreach my $splice (keys %splices) {
      my ($chr,$start,$end,$strand) = split("@",$splice);
      CracTools::SimCT::Utils::printBEDLine($splices_fh,
        { chr     => $chr,
          start   => $start,
          end     => $end, # TODO shoud we add +1 because of half-open bed?
          strand  => $strand,
          name    => join(':',@{$splices{$splice}}),
          score   => scalar @{$splices{$splice}},
        },
      );
    }
  }

  my $info_file   = File::Spec->catfile($simulation_dir,"info.txt");
  my $info_output = CracTools::Output->new(file => $info_file);
  $info_output->printHeaders(
    version => $CracTools::SimCT::VERSION,
    #args    => \@ARGV_copy,
  );
  my %info = (
    nb_reads      => $nb_reads,
    nb_errors     => $nb_errors,
    nb_splices    => scalar keys %splices,
    nb_chimeras   => scalar keys %chimeras,
    nb_mutations  => scalar keys %mutations,
  );
  map { $info_output->printLine($_, $info{$_}) } sort keys %info;

  return $simulation;
}

sub _getReadName {
  my ($read_id, $alignments, $errors_pos) = @_;
  my $read_name = join(":",
    $read_id,
    join(";",
      map {
        join(",",
          $_->chr,
          $_->start,
          $_->strand,
          $_->cigar,
        )
      } @{$alignments},
    ),
  );
  # Encode position list in base 64
  if(@{$errors_pos}) {
    $read_name = join(":", $read_name, CracTools::Utils::encodePosListToBase64(@{$errors_pos}));
  }
  # Truncate read name if greater than 255 char (Maximum length authorized for BAM encoding)
  $read_name = substr $read_name, 0, 255 if length $read_name > 255;
  return $read_name;
}

no Moose;
1;

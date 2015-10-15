package CracTools::SimCT::LiftOver;
# ABSTRACT: Handles coordinates mapping from one reference to an other

use Moose;

use List::Util qw(min max);
use CracTools::Interval::Query;
use CracTools::SimCT::Const;
#use Data::Dumper;

has 'interval_query' => (
  is  => 'ro',
  isa => 'CracTools::Interval::Query',
  default => sub {return CracTools::Interval::Query->new()},
);

# Deletion : offset eq "undef"
sub addInterval {
  my ($self, $chr, $start, $end, $offset, $chr_dest, $reverse) = @_;

  $reverse  = 0    if !defined $reverse;
  $chr_dest = $chr if !defined $chr_dest;

  my $shift_object = {
    offset  => $offset,
    chr     => $chr_dest,
    reverse => $reverse,
  };

  $self->interval_query->addInterval($chr,$start,$end,1,$shift_object);
}

sub shiftInterval {
  my ($self,$chr,$start,$end,$strand) = @_;

  # By default, the strand is positive
  $strand = '+' unless defined $strand;

  my @shifted_intervals;

  # Retrieve intervals from the interval query
  my ($intervals,$offsets) = $self->interval_query->fetchByRegion($chr,$start,$end);

  # Now we sort the intervals indexes
  my @sorted_intervals = sort { $intervals->[$a]->{start} <=> $intervals->[$b]->{start} } (0..@{$intervals}-1);

  foreach my $i (@sorted_intervals) {
    my ($shifted_start, $shifted_end, $ref_start, $ref_end);
    my $offset      = $offsets->[$i]->{offset};
    my $chr_dest    = $offsets->[$i]->{chr};
    my $reverse     = $offsets->[$i]->{reverse};
    my $new_strand  = $strand;

    # Set ref positions
    $ref_start     = max($intervals->[$i]->{start}, $start);
    $ref_end       = min($intervals->[$i]->{end}, $end);

    # If we have an offset, we simply shift coordinates
    if(!$reverse) {
      $shifted_start = $ref_start + $offset;
      $shifted_end   = $ref_end + $offset;
    } else {
      # We have a reverse offset, then all calculous are done
      # the opposite way
      $shifted_end   = $offset - $ref_start;
      $shifted_start = $offset - $ref_end;
      $new_strand    = $new_strand eq '+' ? '-' : '+';
    }

    my $new_interval = {
      start   => $shifted_start,
      end     => $shifted_end,
      ref_start => $ref_start,
      ref_end   => $ref_end,
      chr     => $chr_dest,
      strand  => $new_strand,
      reverse => $reverse,
    };

    push @shifted_intervals, $new_interval;
  }

  return \@shifted_intervals;
}

sub shiftAnnotation {
  my $self = shift;
  my @shifted_intervals = @{$self->shiftInterval(@_)};
  if(@shifted_intervals > 0) {
    if($shifted_intervals[0]->{chr} ne $shifted_intervals[-1]->{chr}) {
      die "Cannot shift annotation overlapping two different chromosomes";
    } 
    if($shifted_intervals[0]->{strand} ne $shifted_intervals[-1]->{strand}) {
      die "Cannot shift annotation overlapping two different strand";
    }
    return {
      chr     => $shifted_intervals[0]->{chr},
      start   => $shifted_intervals[0]->{start},
      end     => $shifted_intervals[-1]->{end},
      strand  => $shifted_intervals[0]->{strand}
    };
  } else {
    return undef;
  }
}

# Given a genomic interval over the simulated genome, 
# returns a set of alignments {start, end, strand, cigar} over the
# reference genome
sub getAlignments {
  my $self = shift;
  my @shifted_intervals = @{$self->shiftInterval(@_)};
  
  # If there is no shifted intervals, we return an empty array
  return () if !@shifted_intervals;

  my $current_interval = shift @shifted_intervals;
  my $match_length = $current_interval->{end} - $current_interval->{start} + 1;
  $current_interval->{cigar} .= $match_length . "M";

  my @alignments = ($current_interval);
  foreach my $shifted_interval (@shifted_intervals) {
    # If this block is collinear with the bed_line genomic coordinate
    # we only update the block_start according to the line start
    if($current_interval->{chr} eq $shifted_interval->{chr} &&
      $current_interval->{strand} eq $shifted_interval->{strand} &&
      $current_interval->{end} < $shifted_interval->{start}) {
      # We have a deletion in the genome we update the cigar with an insertion
      if($shifted_interval->{start} > ($current_interval->{end} + 1)) {
        my $del_length = $shifted_interval->{start} - $current_interval->{end} - 1;
        $current_interval->{cigar} .= $del_length . "D";
      }
      # We have an insertion in the genome we update the cigar with a deletion
      if($shifted_interval->{ref_start} > ($current_interval->{ref_end}+1)) {
        my $ins_length = $shifted_interval->{ref_start} - $current_interval->{ref_end} - 1;
        $current_interval->{cigar} .= $ins_length . "I";
      }
      # Update ref start and end 
      #$current_interval->{ref_start}  = $shifted_interval->{ref_start};
      $current_interval->{ref_end}    = $shifted_interval->{ref_end};
      $current_interval->{end}        = $shifted_interval->{end};
    # Otherwise it is a chimeric alignment
    } else {
      $current_interval = $shifted_interval;
      push @alignments, $current_interval;
    }
    my $match_length = $shifted_interval->{end} - $shifted_interval->{start} + 1;
    $current_interval->{cigar} .= $match_length . "M";
  }
  return @alignments;
}


sub getSplicedAlignments {
  my $self = shift;
  my @intervals = @_;
  my @alignments;
  # Loop over splices
  foreach my $interval (@intervals) {
    # First we get shifted intervals
    my @block_alignments = $self->getAlignments(@{$interval});

    my $prev_alignment = $alignments[$#alignments];
    foreach my $curr_alignment (@block_alignments) {
      # If it is not the first alignment we look for a splice
      # or a chimeric alignement
      if(defined $prev_alignment) {
        # Check if we have a splice alignment
        if($prev_alignment->{chr} eq $curr_alignment->{chr} &&
          $prev_alignment->{strand} eq $curr_alignment->{strand} &&
          $prev_alignment->{reverse} == $curr_alignment->{reverse} &&
          ((!$curr_alignment->{reverse} && 
            $prev_alignment->{end} < $curr_alignment->{start}) ||
          ($curr_alignment->{reverse} &&
            $prev_alignment->{start} > $curr_alignment->{end}))) {
          # Regular spliced alignment
          if(!$curr_alignment->{reverse}) {
            # We can merge the two alignments
            my $splice_length  = $curr_alignment->{start} - $prev_alignment->{end} - 1;
            $prev_alignment->{cigar} .= $splice_length . "N" . $curr_alignment->{cigar};
            $prev_alignment->{end} = $curr_alignment->{end};
          # Reverse spliced alignment
          } else {
            # We can merge the two alignments
            my $splice_length  = $prev_alignment->{start} - $curr_alignment->{end} - 1;
            $prev_alignment->{cigar} = $curr_alignment->{cigar} . $splice_length . "N" . $prev_alignment->{cigar};
            $prev_alignment->{start} = $curr_alignment->{start};
          }
          # Update ref_end pos
          $prev_alignment->{ref_end} = $curr_alignment->{ref_end};
          # Skip this alignement and go to the next one
          next;
        # Otherwise it is a chimeric alignment
        }
      }
      $prev_alignment = $curr_alignment;
      push @alignments, $curr_alignment;
    }
  }
  # Prepend softclip if needed
  my $first_alignment = $alignments[0];
  my $first_interval  = $intervals[0];
  if(defined $first_alignment && $first_alignment->{ref_start} > $first_interval->[1]) {
    my $softclip_length = $first_alignment->{ref_start} - $first_interval->[1];
    $first_alignment->{cigar} = $softclip_length."S".$first_alignment->{cigar};
  }
  # Append softclip if needed
  my $last_alignment = $alignments[$#alignments];
  my $last_interval  = $intervals[$#intervals];
  if(defined $last_alignment && $last_alignment->{ref_end} < $last_interval->[2] - 1) {
    my $softclip_length = $last_interval->[2] - $last_alignment->{ref_end};
    $last_alignment->{cigar} .= $softclip_length."S";
  }
  return @alignments;
}

1;

__END__

=head1 DESCRIPTION

'CracTools::SimCT::LiftOver' is similar to the 'NCBI' LiftOver program by
providing a simple interface for the mapping of genomic coordinates from one
reference to an other.

First, the mapping rule are set by using the 'addInterval' method.  Then the
'shiftInterval' method can be used to realise the map a genomic region to the
other reference.

=head1 METHODS

=head2 new

Create a new 'CracTools::SimCT::LiftOver' object

=head2 addInterval($chr, $start, $end, $offset, $chr_dest, $reverse)

  Arg [1] : 'Str' - Chromosome name of the interval in the old reference
  Arg [2] : 'Int' - Start position of the interval (old ref)
  arg [3] : 'int' - End position of the interval (old ref)
  Arg [4] : 'Int' - Offset to liftover a position from to old reference to the
                    new reference
  Arg [5] : 'Str'     - (Optional) The chromosome name is different on the new
                        reference
  Arg [6] : 'Boolean' - (Optional) The interval is reversed on the new reference

Add a new genomic interval that is defined on the old reference with offset
values to shift this interval to the new reference.

=head2 shiftInterval($chr,$start,$end,$strand) => ArrayRef[HashRef{'Intervals'}]

  Arg [1] : 'Str'    - Chromosome name
  Arg [2] : 'Int'    - Start position of the interval (old ref)
  arg [3] : 'int'    - End position of the interval (old ref)
  Arg [4] : 'Strand' - Strand of the interval

Shift the given interval from the new reference to the old reference and return
a collection of genomic interval. Such an genomic interval is :

  $interval = { chr => '', start => '', end => '', strand => '' }

=head2 shiftAnnotation($chr,$start,$end,$strand)

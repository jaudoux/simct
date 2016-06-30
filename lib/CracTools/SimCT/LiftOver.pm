package CracTools::SimCT::LiftOver;
# ABSTRACT: Handles coordinates mapping from one reference to an other

use Moose;

use List::Util qw(min max);
use CracTools::Interval::Query;
use CracTools::SimCT::Const;
use CracTools::SimCT::GenomicInterval;
use CracTools::SimCT::LiftOver::ShiftedInterval;
use CracTools::SimCT::Alignment;
use CracTools::SimCT::Alignment::CigarElement;
use Data::Dumper;

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
  my ($self,$interval) = @_;

  my @shifted_intervals;

  # Retrieve intervals from the interval query
  my ($intervals,$offsets) = $self->interval_query->fetchByRegion(
    $interval->chr,
    $interval->start,
    $interval->end
  );

  # Now we sort the intervals indexes
  my @sorted_intervals = sort { $intervals->[$a]->{start} <=> $intervals->[$b]->{start} } (0..@{$intervals}-1);

  foreach my $i (@sorted_intervals) {
    my ($shifted_start, $shifted_end, $ref_start, $ref_end);
    my $offset      = $offsets->[$i]->{offset};
    my $chr_dest    = $offsets->[$i]->{chr};
    my $reverse     = $offsets->[$i]->{reverse};
    my $new_strand  = $interval->strand;

    # Set ref positions
    $ref_start     = max($intervals->[$i]->{start}, $interval->start);
    $ref_end       = min($intervals->[$i]->{end}, $interval->end);

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

    my $new_interval = CracTools::SimCT::LiftOver::ShiftedInterval->new(
      chr     => $chr_dest,
      start   => $shifted_start,
      end     => $shifted_end,
      strand  => $new_strand,
      reference_interval => CracTools::SimCT::GenomicInterval->new(
        chr     => $interval->chr,
        start   => $ref_start,
        end     => $ref_end,
        strand  => $interval->strand,
      ),
    );

    push @shifted_intervals, $new_interval;
  }

  return \@shifted_intervals;
}

sub shiftAnnotation {
  my $self = shift;
  my @shifted_intervals = @{$self->shiftInterval(@_)};
  if(@shifted_intervals > 0) {
    if($shifted_intervals[0]->chr ne $shifted_intervals[-1]->chr) {
      die "Cannot shift annotation overlapping two different chromosomes";
    }
    if($shifted_intervals[0]->strand ne $shifted_intervals[-1]->strand) {
      die "Cannot shift annotation overlapping two different strand";
    }
    return CracTools::SimCT::GenomicInterval->new(
      chr     => $shifted_intervals[0]->chr,
      start   => $shifted_intervals[0]->start,
      end     => $shifted_intervals[-1]->end,
      strand  => $shifted_intervals[0]->strand
    );
  } else {
    return undef;
  }
}

# Given a genomic interval over the simulated genome,
# returns a set of alignments {start, end, strand, cigar} over the
# reference genome
sub getAlignments {
  my $self = shift;
  my $original_interval = shift;

  my @shifted_intervals = @{$self->shiftInterval($original_interval)};

  # If there is no shifted intervals, we return an empty array
  return () if !@shifted_intervals;

  my $first_interval = shift @shifted_intervals;

  # Create a first alignement object
  my $query_mapping_start = $first_interval->reference_interval->start - $original_interval->start;
  my $current_alignment = CracTools::SimCT::Alignment->new(
    chr     => $first_interval->chr,
    start   => $first_interval->start,
    strand  => $first_interval->strand,
    query_length  => $original_interval->length,
    query_mapping_start => $query_mapping_start,
  );

  # Append the first interval to the alignment
  $current_alignment->appendCigarElement(
    CracTools::SimCT::Alignment::CigarElement->new(
      op => 'M',
      nb => $first_interval->length,
    ),
  );

  my @alignments = ($current_alignment);

  # TODO change prev_interval by "current_alignment"
  my $prev_interval = $first_interval;
  foreach my $shifted_interval (@shifted_intervals) {
    # If this block is collinear with the bed_line genomic coordinate
    # we only update the block_start according to the line start
    if($prev_interval->chr eq $shifted_interval->chr &&
      $prev_interval->strand eq $shifted_interval->strand &&
      $prev_interval->end < $shifted_interval->start) {
      # We have a deletion in the genome we update the cigar with an insertion
      if($shifted_interval->start > ($prev_interval->end + 1)) {
        $current_alignment->appendCigarElement(
          CracTools::SimCT::Alignment::CigarElement->new(
            op => 'D',
            nb => $shifted_interval->start - $prev_interval->end - 1,
          ),
        );
      }
      # We have an insertion in the genome we update the cigar with a deletion
      if($shifted_interval->reference_interval->start > ($prev_interval->reference_interval->end+1)) {
        $current_alignment->appendCigarElement(
          CracTools::SimCT::Alignment::CigarElement->new(
            op => 'I',
            nb => $shifted_interval->reference_interval->start - $prev_interval->reference_interval->end - 1,
          ),
        );
      }
    # Otherwise it is a chimeric alignment
    } else {
      # Create a new alignement object
      $query_mapping_start = $shifted_interval->reference_interval->start - $original_interval->start;
      $current_alignment = CracTools::SimCT::Alignment->new(
        chr     => $shifted_interval->chr,
        start   => $shifted_interval->start,
        strand  => $shifted_interval->strand,
        query_length  => $original_interval->length,
        query_mapping_start => $query_mapping_start,
      );
      push @alignments, $current_alignment;
    }
    $current_alignment->appendCigarElement(
      CracTools::SimCT::Alignment::CigarElement->new(
        op => 'M',
        nb => $shifted_interval->length,
      ),
    );
    $prev_interval = $shifted_interval;
  }
  return @alignments;
}


sub getSplicedAlignments {
  my $self = shift;
  my @intervals = @_;
  my @alignments;

  # Compute the whole length of the query
  my $query_length = 0;
  map { $query_length += $_->length } @intervals;

  # Loop over splices
  foreach my $interval (@intervals) {
    # First we get shifted intervals
    my @block_alignments = $self->getAlignments($interval);

    my $prev_alignment = $alignments[$#alignments];
    foreach my $curr_alignment (@block_alignments) {
      # If it is not the first alignment we look for a splice
      # alignement
      if(defined $prev_alignment &&
        $prev_alignment->chr eq $curr_alignment->chr &&
        $prev_alignment->strand eq $curr_alignment->strand &&
        $prev_alignment->end < $curr_alignment->start) {

        # We merge two alignments, we need to transform the softclips into insertions
        if($prev_alignment->right_softclip > 0) {
          $prev_alignment->appendCigarElement(CracTools::SimCT::Alignment::CigarElement->new(
            op => 'I',
            nb => $prev_alignment->right_softclip,
          ));
        }
        # Regular spliced alignment
        # We can merge the two alignments
        my $splice_length  = $curr_alignment->start - $prev_alignment->end - 1;
        $prev_alignment->appendCigarElement(CracTools::SimCT::Alignment::CigarElement->new(
          op => 'N',
          nb => $splice_length,
        ));
        foreach my $cigel ($curr_alignment->allCigarElements) {
          $prev_alignment->appendCigarElement($cigel);
        }
        # Now we update the query length
        $prev_alignment->query_length($query_length);
        #next;
      # Otherwise this is a chimeric alignment and we push this alignement
      } else {
        if(defined $prev_alignment) {
          $prev_alignment->query_length($query_length);
          $curr_alignment->query_mapping_start($prev_alignment->query_mapping_end + 1);
        }
        $prev_alignment = $curr_alignment;
        push @alignments, $curr_alignment;
      }
    }
  }
  return @alignments;
}

no Moose;
__PACKAGE__->meta->make_immutable;

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

  Arg [1] : 'CracTools::SimCT::GenomicInterval'  - Genomic interval to shift

Shift the given interval from the new reference to the old reference and return
a list of CracTools::SimCT::GenomicInterval.

=head2 shiftAnnotation($chr,$start,$end,$strand)

package CracTools::SimCT::LiftOver;

use Moose;

use List::Util qw(min max);
use CracTools::Interval::Query;

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

# strand => '1|-1'
sub shiftInterval {
  my ($self,$chr,$start,$end,$strand) = @_;

  # By default, the strand is positive
  $strand = '+' unless defined $strand;

  my @shifted_intervals;

  # Retrieve intervals from the interval query
  my ($intervals,$offsets) = $self->interval_query->fetchByRegion($chr,$start,$end);

  # Now we sort the intervals
  my @sorted_intervals = sort { $intervals->[$a]->{start} <=> $intervals->[$b]->{start} } (0..@{$intervals}-1);

  foreach my $i (@sorted_intervals) {
    my ($shifted_start, $shifted_end);
    my $offset      = $offsets->[$i]->{offset};
    my $chr_dest    = $offsets->[$i]->{chr};
    my $reverse     = $offsets->[$i]->{reverse};
    my $new_strand  = $strand;

    # If we have an offset, we simply shift coordinates
    if(!$reverse) {
      $shifted_start = $start + $offset;
      $shifted_end   = min($intervals->[$i]->{end} + $offset, $end + $offset);
    } else {
      # We have a reverse offset, then all calculous are done
      # the opposite way
      $shifted_end   = $offset - $start;
      $shifted_start = min($offset - $intervals->[$i]->{start}, $offset - $end);
      $new_strand    = $new_strand eq '+' ? '-' : '+';
    }

    my $new_interval = {
      start   => $shifted_start,
      end     => $shifted_end,
      chr     => $chr_dest,
      strand  => $new_strand,
    };

    $start += $new_interval->{end} - $new_interval->{start} + 1;
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

1;

package CracTools::SimCT::Fusion::FusedExon;

use Moose::Role;

use CracTools::Utils;
use CracTools::SimCT::Annotations::Exon;

has 'exon' => (
  is       => 'ro',
  isa      => 'CracTools::SimCT::Annotations::Exon',
  required => 1,
);

has 'fused_sequence' => (
  is       => 'rw',
  isa      => 'DNA',
);

requires 'setFusedSequence';

requires 'allFusedExons';

# Fast accessor for chr
sub chr {
  my $self = shift;
  return $self->gene->chr;
}

# Fast accessor for strand
sub strand {
  my $self = shift;
  return $self->gene->strand;
}

# Fast accessor for gene
sub gene {
  my $self = shift;
  return $self->exon->gene;
}

# TODO this method should work even if 'fused_sequence' is not
# defined yet using the exon positions
sub length {
  my $self = shift;
  return length $self->fused_sequence;
}

package CracTools::SimCT::Fusion::FusedExon::5prim;

use Moose;
use Scalar::Util 'refaddr';

with 'CracTools::SimCT::Fusion::FusedExon';

sub setFusedSequence {
  my $self = shift;
  my $chr_seq_ref = shift;
  my $fused_seq;
  if($self->strand eq '+') {
    $fused_seq = substr ${$chr_seq_ref}, $self->gene->start, $self->exon->end - $self->gene->start + 1;
  } else {
    $fused_seq = substr ${$chr_seq_ref}, $self->exon->start, $self->gene->end - $self->exon->start + 1;
    $fused_seq = CracTools::Utils::reverseComplement($fused_seq);
  }
  $self->fused_sequence($fused_seq);
}

# Return all exon to be fused, the last one beeing the fused exon
sub allFusedExons {
  my $self = shift;
  my @fused_exons;
  if($self->strand eq '+') {
    my @exons = $self->gene->sortedExons;
    foreach my $exon (@exons) {
      if(refaddr($exon) == refaddr($self->exon)) {
        push @fused_exons, $exon; 
        last;
      } elsif($exon->end < $self->exon->end) {
        push @fused_exons, $exon; 
      }
    }
  } else {
    my @exons = sort { $b->end <=> $a->end } $self->gene->allExons;
    foreach my $exon (@exons) {
      if(refaddr($exon) == refaddr($self->exon)) {
        push @fused_exons, $exon; 
        last; 
      } elsif($exon->start > $self->exon->start) {
        push @fused_exons, $exon; 
      }
    }
  }
  return @fused_exons;
}

package CracTools::SimCT::Fusion::FusedExon::3prim;

use Moose;
use Scalar::Util 'refaddr';

with 'CracTools::SimCT::Fusion::FusedExon';

sub setFusedSequence {
  my $self = shift;
  my $chr_seq_ref = shift;
  my $fused_seq;
  if($self->strand eq '+') {
    $fused_seq = substr ${$chr_seq_ref}, $self->exon->start, $self->gene->end - $self->exon->start + 1;
  } else {
    $fused_seq = substr ${$chr_seq_ref}, $self->gene->start, $self->exon->end - $self->gene->start + 1;
    $fused_seq = CracTools::Utils::reverseComplement($fused_seq);
  }
  $self->fused_sequence($fused_seq);
}

# Return all exon to be fused, the first one being the fused exon
sub allFusedExons {
  my $self = shift;
  my @fused_exons;
  if($self->strand eq '+') {
    my @exons = $self->gene->sortedExons;
    foreach my $exon (@exons) {
      if(refaddr($exon) == refaddr($self->exon)) {
        @fused_exons = ();
        push @fused_exons, $exon; 
      } elsif($exon->start > $self->exon->start) {
        push @fused_exons, $exon; 
      }
    }
  } else {
    my @exons = sort { $b->end <=> $a->end } $self->gene->allExons;
    foreach my $exon (@exons) {
      if(refaddr($exon) == refaddr($self->exon)) {
        @fused_exons = ();
        push @fused_exons, $exon; 
      } elsif($exon->end < $self->exon->end) {
        push @fused_exons, $exon; 
      }
    }
  }
  return @fused_exons;
}

1;

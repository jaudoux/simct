package CracTools::SimCT::Fusion::FusedExon;
# ABSTRACT: A fused exon implicated in a 'CracTools::SimCT::Fusion'

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

__END__

=head1 DESCRIPTION

'CracTools::SimCT::Fusion::FusedExon' holds an exon implicated in a gene fusion.

This base class has two implementations :

=over 2

=item - L<CracTools::SimCT::Fusion::FusedExon::5prim>

Which describes a fused exon that implicates the 5prim extremity of its gene

=item - L<CracTools::SimCT::Fusion::FusedExon::3prim>

Which describes a fused exon that implicates the 3prim extremity of its gene

=back

=head1 ACCESSORS

=head2 exon => 'CracTools::SimCT::Annotations::Exon'

Getter for the fused exon

=head2 fused_sequence => 'DNA'

Getter for the fused sequence

=head2 length

Length of the fused sequence

=head2 gene => 'CracTools::SimCT::Annotations::Gene'

Fast accessor for the exon's gene (eq. $self->exon->gene)

=head2 chr => 'Str'

Fast accessor for the exon's chromosome

=head2 strand => 'Strand'

Fast accessor for the exon's strand

=head1 METHODS

=head2 new

  Arg [exon] : 'CracTools::SimCT::Annotations::Exon' - the fused exon

Create a new 'CracTools::SimCT::Fusion::FusedExon' object

=head2 setFusedSequence($chr_ref)

Given a reference of the chromosome sequence, set the 'fused_sequence'
attributes.

For the 'CracTools::SimCT::Fusion::FusedExon::5prim' it correspond to the
5prim start of the gene until the end of the fused exon.

For the 'CracTools::SimCT::Fusion::FusedExon::3prim' it correspond to the
start of the fused exon until the 3prim end of the gene.

=head2 allFusedExons => Array('CracTools::SimCT::Annotations::Exons')

Return an array of all implicated exons in the fusion for this fused gene.

In the case of the 'CracTools::SimCT::Fusion::FusedExon::5prim' implementation,
the fused exon, is stored in the last position of the array (ie. accessible with
a 'pop').

In the case of the 'CracTools::SimCT::Fusion::FusedExon::3prim' implementation,
the fused exon, is stored in the first position of the array (ie. accessible with
a 'shift').


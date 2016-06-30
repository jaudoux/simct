package CracTools::SimCT::Fusion;
# ABSTRACT: A fusion gene formed by two fused exons

use Moose;
use Scalar::Util 'refaddr';

use CracTools::SimCT::Fusion::FusedExon;
use CracTools::SimCT::Const;

#with ('CracTools::SimCT::GeneticVariant');

has 'fused_exon_5prim' => (
  is        => 'ro',
  isa       => 'CracTools::SimCT::Fusion::FusedExon::5prim',
  required  => 1,
);

has 'fused_exon_3prim' => (
  is        => 'ro',
  isa       => 'CracTools::SimCT::Fusion::FusedExon::3prim',
  required  => 1,
);

has 'chr_fusion_pos' => (
  is => 'rw',
  isa => 'Int',
  default => 0,
);

has 'chr_fusion' => (
  is => 'rw',
  isa => 'Str',
  default => "$CracTools::SimCT::Const::CHR_FUSIONS",
);

sub fusion_sequence {
  my $self = shift;
  return $self->fused_exon_5prim->fused_sequence.$self->fused_exon_3prim->fused_sequence;
}

sub isFusionSequenceDefined {
  my $self = shift;
  if(defined $self->fused_exon_5prim->fused_sequence &&
     defined $self->fused_exon_3prim->fused_sequence) {
    return 1;
  } else {
    return 0;
  }
}

sub getFusionGene($$$) {
  my $self = shift;
  my ($fusion_id,$chr_fusion,$chr_fusion_pos) = @_;

  $self->chr_fusion_pos($chr_fusion_pos);
  $self->chr_fusion($chr_fusion);

  # Create a new fusion gene object
  my $fusion_gene = CracTools::SimCT::Annotations::Gene->new(
    chr     => $chr_fusion,
    id      => $fusion_id,
    strand  => '+',
  );

  # Create the fusion gene exons
  my @fusion_exons;
  foreach my $fused_exon (($self->fused_exon_5prim,$self->fused_exon_3prim)) {
    my @exons = $fused_exon->allFusedExons;
    # If this is the 3prim fused exon, we fused it with the 5prim exon
    if(refaddr($fused_exon) == refaddr($self->fused_exon_3prim)) {
      my $fused_exon_3prim  = shift @exons;
      my $last_fusion_exon  = $fusion_exons[$#fusion_exons];
      $chr_fusion_pos = $last_fusion_exon->end + 1;
      # Fuse the 3prim exon with the 5prim exon into a uniq annotation
      $last_fusion_exon->end($last_fusion_exon->end + $fused_exon_3prim->length);
    }
    # We loop over the exons and create the corresponding annotations
    foreach my $exon (@exons) {
      my $exon_start;
      if(refaddr($fused_exon) == refaddr($self->fused_exon_5prim)) {
        $exon_start = $exon->strand eq '+'?
                       $chr_fusion_pos + ($exon->start - $exon->gene->start) :
                       $chr_fusion_pos + ($exon->gene->end - $exon->end);
      } else {
        $exon_start = $exon->strand eq '+'?
                       $chr_fusion_pos + ($exon->start - $fused_exon->exon->start) :
                       $chr_fusion_pos + ($fused_exon->exon->end - $exon->end);
      }
      push @fusion_exons, CracTools::SimCT::Annotations::Exon->new(
        gene  => $fusion_gene,
        start => $exon_start,
        end   => $exon_start + $exon->length - 1,
        transcripts  => [$fusion_id.".1"],
      );
    }
  }

  return $fusion_gene;
}

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 DESCRIPTION

A fusion is a chimeric gene formed by fusing the exons of two genes (not
necessarly different). One of the implicated exon is said 5prim and the other 3prim.
The new fused gene will be formed by fusing the upstream exon's of the 5prim gene
and the downstream exons of the 3prim gene. The two fused exons will be merge
into a new exon.

=head1 ACCESSORS

=head2 fused_exon_5prim => 'CracTools::SimCT::Fusion::FusedExon::5prim'

Getter for the 5prim fused exon

=head2 fused_exon_3prim => 'CracTools::SimCT::Fusion::FusedExon::3prim'

Getter for the 3prim fused exon

=head2 fused_sequence => 'DNA'

Getter for the full gene fusion sequence wich is the merge of the two
fused gene's fusion sequences.

=head1 METHODS

=head2 new

  Arg [fused_exon_5prim] : 'CracTools::SimCT::Fusion::FusedExon::5prim' - the
                            5prim fused exon
  Arg [fused_exon_3prim] : 'CracTools::SimCT::Fusion::FusedExon::3prim' - the
                            3prim fused exon

Create a new 'CracTools::SimCT::Fusion' object

=head2 isFusionSequenceDefined => 1|0

Return true if the fusion sequence is completely defined (ie. both 5prim and
3prim sequences have been defined)

=head2 getFusionGene => 'CracTools::SimCT::Annotations::Gene'

  Arg [1] : 'Str' - Fusion id used to name the fusion gene
  Arg [2] : 'Str' - Chromosome name
  Arg [3] : 'Int' - Chromosome position where the fusion start

Return a new 'CracTools::SimCT::Annotations::Gene' object that correspond
to the fusion gene.

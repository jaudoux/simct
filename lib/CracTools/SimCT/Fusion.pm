package CracTools::SimCT::Fusion;

use Moose;
use CracTools::SimCT::Fusion::FusedExon;
use CracTools::SimCT::Const;

with ('CracTools::SimCT::GeneticVariant');

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

has 'fusion_id' => (
  is        => 'rw',
  isa       => 'Str',
  default   => 'fusion',
);

has 'chr_fusion' => (
  is        => 'rw',
  isa       => 'Str',
  default   => "$CracTools::SimCT::Const::CHR_FUSIONS",
);

has 'chr_fusion_pos' => (
  is        => 'rw',
  isa       => 'Int',
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

sub getFusionGene {
  my $self = shift;

  my $chr_fusion_pos = $self->chr_fusion_pos;
  my $fused_exon;

  die "Cannot create fusion gene because chr_fusion_pos is not defined"
  unless defined $self->chr_fusion_pos;

  my $fusion_gene = CracTools::SimCT::Annotations::Gene->new(
    chr     => $self->chr_fusion,
    id      => $self->fusion_id,
    strand  => '+',
  );
  

  my @fused_exons_5prim = $self->fused_exon_5prim->allFusedExons;

  foreach my $exon (@fused_exons_5prim) {
    my $exon_start = $exon->strand eq '+'? 
                     $chr_fusion_pos + ($exon->start - $exon->gene->start) : 
                     $chr_fusion_pos + ($exon->gene->end - $exon->end);
    my $new_exon  = CracTools::SimCT::Annotations::Exon->new(
      gene  => $fusion_gene,
      start => $exon_start,
      end   => $exon_start + $exon->length - 1,
      transcripts  => [$self->fusion_id.".1"],
    );
    $fused_exon = $new_exon;
  }

  $chr_fusion_pos = $fused_exon->end + 1;

  my @fused_exons_3prim = $self->fused_exon_3prim->allFusedExons;
  my $fused_exon_3prim  = shift @fused_exons_3prim;

  # Fuse the 3prim exon with the 5prim exon into a uniq annotation
  $fused_exon->end($fused_exon->end + $fused_exon_3prim->length);

  foreach my $exon (@fused_exons_3prim) {
    my $exon_start = $exon->strand eq '+'? 
                     $chr_fusion_pos + ($exon->start - $exon->gene->start) : 
                     $chr_fusion_pos + ($exon->gene->end - $exon->end);
    my $new_exon  = CracTools::SimCT::Annotations::Exon->new(
      gene  => $fusion_gene,
      start => $exon_start,
      end   => $exon_start + $exon->length - 1,
      transcripts  => [$self->fusion_id.".1"],
    );
  }

  return $fusion_gene;
}

1;

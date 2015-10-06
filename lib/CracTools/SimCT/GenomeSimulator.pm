package CracTools::SimCT::GenomeSimulator;
# ABSTRACT: Simulate a mutated genome by introducing indel, snv and fusion genes in a reference genome

use Moose;

use CracTools::GenomeMask;

use CracTools::SimCT::Utils;
use CracTools::SimCT::LiftOver;
use CracTools::SimCT::Annotations;
use CracTools::SimCT::SimulatedGenome;

use Carp;
use File::Spec;

has 'genome' => (
  is => 'ro',
  isa => 'CracTools::SimCT::Genome',
  required => 1,
);

# This genome mask keeps track of mutation position
# that have already been added in order to avoid
# conflict
has 'genome_mask' => (
  is => 'ro',
  isa => 'CracTools::GenomeMask',
  lazy => 1,
  default => sub {
    my $self = shift;
    my $genome_mask = CracTools::GenomeMask->new(
      genome => $self->genome->references_length,
    );
    return $genome_mask;
  },
);

# Hold the mutations
has 'mutations' => (
  traits => ['Array'],
  is => 'rw',
  isa => 'ArrayRef[CracTools::SimCT::Mutation]',
  default => sub { [] },
  handles => {
    allMutations => 'elements',
  },
);

has 'fusions' => (
  traits => ['Array'],
  is => 'rw',
  isa => 'ArrayRef[CracTools::SimCT::Fusion]',
  default => sub { [] },
  handles => {
    allFusions  => 'elements',
    addFusion   => 'push',
  },
);

# Return true is the mutation has been properly added
sub addMutation($) {
  my $self = shift;
  my $mut = shift;
  my $ref_length = $self->genome->getReferenceLength($mut->chr);
  if(!defined $ref_length) {
    carp "Mutation chromosome is not available in the reference genome";
  } elsif ($mut->end >= $self->genome->getReferenceLength($mut->chr)) {
    carp "Mutation is outside the scope of the reference";
  } elsif($self->genome_mask->getNbBitsSetInRegion($mut->chr,$mut->pos,$mut->end+1) == 0) {
    push @{$self->mutations}, $mut;
    $self->genome_mask->setRegion($mut->chr,$mut->pos,$mut->end+1);
    return 1;
  }
  return 0;
}

sub sortedMutations {
  my $self = shift;
  my $chr  = shift;
  if(defined $chr) {
    return sort {$a->pos <=> $b->pos} grep {$_->chr eq $chr} $self->allMutations;
  } else {
    return sort {$a->chr cmp $b->chr && $a->pos <=> $b->pos} $self->allMutations;
  }
}

sub generateGenome {
  my $self = shift;
  my %args = @_;
  
  my $annotations = $args{annotations};
  my $genome_dir  = $args{genome_dir};

  my $gtf_file        = File::Spec->catfile($genome_dir,"annotations.gtf");
  my $gtf_output_fh   = CracTools::Utils::getWritingFileHandle($gtf_file);

  my $annot_lifter = CracTools::SimCT::LiftOver->new();
  
  # We generate the simulated genome for each chr
  foreach my $chr ($self->genome->sortedReferences) {
    my $fasta_output    = File::Spec->catfile($genome_dir,"$chr.fa");
    my $fasta_output_fh = CracTools::Utils::getWritingFileHandle($fasta_output);
    my $chr_seq         = $self->genome->getReferenceSeq($chr);
    my $chr_length      = $self->genome->getReferenceLength($chr);
  
    # We update the fusion genes that involves this chromosome and
    # store the corresponding sequences
    foreach my $fusion ($self->allFusions) {
      if($fusion->fused_exon_5prim->chr eq $chr) {
        $fusion->fused_exon_5prim->setFusedSequence(\$chr_seq);
      }
      if($fusion->fused_exon_3prim->chr eq $chr) {
        $fusion->fused_exon_3prim->setFusedSequence(\$chr_seq);
      }
    }

    my $remainder       = 0; # What's left in the current FASTA line
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome

    # Write the header of fasta output
    print $fasta_output_fh ">$chr\n";

    foreach my $mut ($self->sortedMutations($chr)) {

      # Add the interval between the previous mutation and the current one
      # to the liftover
      $annot_lifter->addInterval($chr,$index,$mut->pos-1,$offset);

      # We read the sequence before the mutation and print it in
      # the output fasta
      #my $frag    = substr $chr_seq, $index, $mut->pos-$index, "";
      my $frag    = substr $chr_seq, 0, $mut->pos-$index, "";
      $remainder  = CracTools::SimCT::Utils::printFASTA($fasta_output_fh,$frag,$remainder);

      # Give the a reference to the chr sequence for the mutation
      # to load its reference sequence
      my $mut_ref_seq = substr $chr_seq, 0, $mut->referenceLength, "";
      $mut->setReferenceSequence($mut_ref_seq);

      # Print the sequence that correspond to the mutation in
      # the output FASTA
      $remainder = CracTools::SimCT::Utils::printFASTA(
        $fasta_output_fh,
        $mut->mutation_sequence,
        $remainder,
      );

      # Update the offset and the index
      $offset += ($mut->mutationLength - $mut->referenceLength);
      $index  = $mut->pos + $mut->referenceLength;
    }

    # Now we print what's left of the reference
    my $frag    = substr $chr_seq, 0, $chr_length - $index;
    CracTools::SimCT::Utils::printFASTA($fasta_output_fh,$frag,$remainder);

    # Add the last interval to the liftover
    $annot_lifter->addInterval($chr,$index,$chr_length-1,$offset);


    close($fasta_output_fh);
  }

  # Remove silencious mutations
  my @cleaned_mutations = grep { $_->reference_sequence ne $_->mutation_sequence } $self->allMutations;
  $self->mutations($self->mutations);

  # Remove undefined fusions
  my @cleaned_fusions = grep { $_->isFusionSequenceDefined } $self->allFusions;
  $self->fusions(\@cleaned_fusions);

  # Print liftover annotations
  $annotations->appendGTF($gtf_output_fh,$annot_lifter);

  # Now we print an extra FASTA file with fusions
  my $fasta_output    = File::Spec->catfile($genome_dir,"chr$CracTools::SimCT::Const::CHR_FUSIONS.fa");
  my $fasta_output_fh = CracTools::Utils::getWritingFileHandle($fasta_output);

  # Print FASTA headers
  print $fasta_output_fh ">chr$CracTools::SimCT::Const::CHR_FUSIONS\n";

  my $remainder       = 0;
  my $fusion_id       = 0;
  my $chr_fusion_pos  = 0;

  # Create a new annotation set for fusions
  my $fusion_annotations = CracTools::SimCT::Annotations->new();

  foreach my $fusion ($self->allFusions) {
    # We print the fasta sequence corresponding to the fusion
    $remainder = CracTools::SimCT::Utils::printFASTA($fasta_output_fh,$fusion->fusion_sequence,$remainder);

    # No we update the annotations and print them to the GTF
    # Update the position of the fusion in the "$CracTools::SimCT::Const::CHR_FUSIONS" chr
    # TODO this could be parameters of the 'getFusionGene' method !?
    $fusion->fusion_id("fusion_$fusion_id");
    $fusion->chr_fusion($CracTools::SimCT::Const::CHR_FUSIONS);
    $fusion->chr_fusion_pos($chr_fusion_pos);

    # Get the new genes and add it to the annotations
    my $fusion_gene = $fusion->getFusionGene;
    $fusion_annotations->addGene($fusion_gene);
    $chr_fusion_pos = $fusion_gene->end + 1;
    $fusion_id++;
  }

  # Print fusion annotations
  $fusion_annotations->appendGTF($gtf_output_fh);

  # Close outputs
  close($fasta_output_fh);
  close($gtf_output_fh);
  
  # Return the simulated genome build on the base
  # of the current state of GS
  return CracTools::SimCT::SimulatedGenome->new(
    genome_simulator => $self,
  );
}

1;

__END__

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
  } elsif($self->genome_mask->getNbBitsSetInRegion($mut->chr,$mut->start,$mut->end+1) == 0) {
    push @{$self->mutations}, $mut;
    $self->genome_mask->setRegion($mut->chr,$mut->start,$mut->end+1);
    return 1;
  }
  return 0;
}

sub sortedMutations {
  my $self = shift;
  my $chr  = shift;
  if(defined $chr) {
    return sort {$a->start <=> $b->start} grep {$_->chr eq $chr} $self->allMutations;
  } else {
    return sort {$a->chr cmp $b->chr || $a->start <=> $b->start} $self->allMutations;
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

    my $prev_pos        = 0;
    my $remainder       = 0; # What's left in the current FASTA line
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome

    # Write the header of fasta output
    print $fasta_output_fh ">$chr\n";

    foreach my $mut ($self->sortedMutations($chr)) {

      # We read the sequence before the mutation and print it in
      # the output fasta
      my $frag    = substr $chr_seq, 0, $mut->start - $index, "";
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

      # Update prev_pos and index
      $prev_pos = $index;
      $index    = $mut->start + CracTools::SimCT::Utils::min($mut->referenceLength,$mut->mutationLength);

      # Add the current interval
      $annot_lifter->addInterval($chr,$prev_pos,$index-1,$offset);

      # Update the offset and the index
      $offset += ($mut->mutationLength - $mut->referenceLength);
      $index  = $mut->start + $mut->referenceLength;

    }

    # Now we print what's left of the reference
    my $frag = substr $chr_seq, 0, $chr_length - $index;
    CracTools::SimCT::Utils::printFASTA($fasta_output_fh,$frag,$remainder);

    # Add the last interval to the liftover
    if($index < $chr_length) {
      $annot_lifter->addInterval($chr,$index,$chr_length-1,$offset);
    }

    close($fasta_output_fh);
  }

  # Remove silencious mutations
  my @cleaned_mutations = grep { $_->reference_sequence ne $_->mutation_sequence } $self->allMutations;
  $self->mutations(\@cleaned_mutations);

  # Remove undefined fusions
  my @cleaned_fusions = grep { $_->isFusionSequenceDefined } $self->allFusions;
  $self->fusions(\@cleaned_fusions);

  # Print liftover annotations
  $annotations->appendGTF($gtf_output_fh,$annot_lifter);


  my $fusion_id = 0;

  # Create a new annotation set for fusions
  my $fusion_annotations = CracTools::SimCT::Annotations->new();

  foreach my $fusion ($self->allFusions) {

    my $fusion_name = "fusion_$fusion_id";

    # Now we print an extra FASTA file with fusions
    my $fasta_output    = File::Spec->catfile($genome_dir,"$fusion_name.fa");
    my $fasta_output_fh = CracTools::Utils::getWritingFileHandle($fasta_output);

    # Print FASTA headers
    print $fasta_output_fh ">$fusion_name\n";
    # We print the fasta sequence corresponding to the fusion
    CracTools::SimCT::Utils::printFASTA($fasta_output_fh,$fusion->fusion_sequence,0);

    # Get the new genes and add it to the annotations
    my $fusion_gene = $fusion->getFusionGene(
      $fusion_name,
      $fusion_name,
      0,
    );
    $fusion_annotations->addGene($fusion_gene);

    $fusion_id++;
    close($fasta_output_fh);
  }

  # Print fusion annotations
  $fusion_annotations->appendGTF($gtf_output_fh);

  # Close outputs
  close($gtf_output_fh);

  # Return the simulated genome build on the base
  # of the current state of GS
  return ($gtf_file, CracTools::SimCT::SimulatedGenome->new(
    genome_simulator => $self,
  ));
}

no Moose;
__PACKAGE__->meta->make_immutable;

__END__

=head1 DESCRIPTION

L<CracTools::SimCT::GenomeSimulator> is a fundamental class of SimCT wich will
create a L<CracTools::SimCT::SimulatedGenome> given a reference genome and a set of
mutations and fusions. Multiple simulated genome can be generated from the same
GenomeSimulator by introducing new mutations and/or removing some of them.

=head1 ACCESSORS

=head2 genome => CracTools::SimCT::Genome

Getter for the reference genome

=head2 genome_mask => CracTools::GenomeMask

Getter for the genome mask associated to the reference genome

=head2 mutations => ArrayRef[CracTools::SimCT::Mutation]

Getter for the array references that contains the mutations introduced in the
simulated genome.

=head1 METHODS

=head2 new

  Arg [genome]  : 'CracTools::SimCT::Genome' - the reference genome to be used

Create a new 'CracTools::SimCT::GenomeSimulator' object

=head2 allMutations => Array('CracTools::SimCT::Mutation')

Return all inserted mutations

=head2 sortedMutations  => Array('CracTools::SimCT::Mutation')

Return all inserted mutations sorted by chromosome and position

=head2 addMutation('CracTools::SimCT::Mutation') => 0|1

Add a mutations and return true if the mutation was properly added or false if
it was in conflict with a previously inserted mutation or outside the scope of
the reference genome (ie. chromosome sequence not available or position greater
than the chromosome length).

=head2 allFusions => Array('CracTools::SimCT::Fusion')

Return all inserted fusions

=head2 addFusion('CracTools::SimCT::Fusion')

Add a fusion

=head2 generateGenome => 'CracTools::SimCT::SimulatedGenome'

  Arg [genome_dir]   : 'Path' - the directory where the simulated genome will be
                        saved.
  Arg [annotations]  : 'CracTools::SimCT::Annotations' - A set of annotations
                        over the reference genome that will be liftover the
                        simulated genome and print into
                        "genome_dir/annotations.gtf"

Simulated a genome using the current state (ie. introduced mutations) of the
GenomeSimulator. Genome chromosome sequences are saved in FASTA format in
separated files. If an annotations object is provided they are shifted over the
simulated genome and also saved in the 'GTF' format in the same directory under
'annotations.gtf' filename.

This method also return a 'CracTools::SimCT::SimulatedGenome' new object wich
is able shift chromosomic coordinates from the simulated genome to the reference
genome. 'CracTools::SimCT::SimulatedGenome' also hold a frozen set of mutations
that where used to generate the simulated genome.

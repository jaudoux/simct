package CracTools::SimCT::SimulatedGenome;
# ABSTRACT: Represent a Simulated genome and contains mapping values between the original genome

use Moose;

require 5.004;
use Tie::RefHash;

use CracTools::Utils;
use CracTools::SimCT::Utils;

has 'liftover' => (
  is  => 'ro',
  isa => 'CracTools::SimCT::LiftOver',
  default => sub { CracTools::SimCT::LiftOver->new(); },
);

has 'mutation_query' => (
  is  => 'ro',
  isa => 'CracTools::Interval::Query',
  default => sub { CracTools::Interval::Query->new(); },
);

sub BUILD {
  my $self = shift;
  my $args = shift;
  my $genome_simulator = $args->{genome_simulator};

  # Create a liftover object
  foreach my $chr (sort $genome_simulator->genome->references) {
    my $prev_pos        = 0;
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome
    my $chr_length      = $genome_simulator->genome->getReferenceLength($chr);

    foreach my $mut ($genome_simulator->sortedMutations($chr)) {
      
      $prev_pos = $index;
      $index    = $mut->pos - $offset;

      # Add the current interval
      $self->liftover->addInterval($chr,$prev_pos,$index-1,$offset);

      # Add the mutation into the mutation query
      $self->mutation_query->addInterval($chr,$index,$index,1,$mut);

      # Update offset and index
      $index  += $mut->mutationLength; 
      $offset += ($mut->referenceLength - $mut->mutationLength);
    }
    # Add the last interval
    $self->liftover->addInterval($chr,$index,$chr_length-1,$offset);
  }

  # Add fusions in the liftover query
  foreach my $fusion ($genome_simulator->allFusions) {
    my $start_pos   = $fusion->chr_fusion_pos;
    foreach my $fused_exon (($fusion->fused_exon_5prim,$fusion->fused_exon_3prim)) {
      # Add interval for the 5prim part
      $self->liftover->addInterval(
        $fusion->chr_fusion,
        $start_pos,
        $start_pos + $fused_exon->length - 1,
        $fused_exon->strand eq '+'?
          $fused_exon->gene->start  - $start_pos :
          $fused_exon->gene->end    + $start_pos,
        $fused_exon->chr,
        $fused_exon->strand eq '+'? 0 : 1,
      );
      # Update the start pos for 3prim fused exon
      $start_pos    += $fused_exon->length;
    }
  }
};


sub shiftAlignements {
  my $self      = shift;
  my %args      = @_;

  my $bed_file          = $args{bed_file};
  my $output_file       = $args{output_file};
  my $vcf_file          = $args{vcf_file};
  my $chimera_file      = $args{chimera_file};

  my $output_fh         = CracTools::Utils::getWritingFileHandle($output_file);
  my $vcf_fh            = CracTools::Utils::getWritingFileHandle($vcf_file);

  # TODO hand chimera output file

  # Hold mutations and count, we use 'Tie::RefHash' to use mutations (hash)
  # as hash keys
  tie my %overlapped_mutations, 'Tie::RefHash';

  # Now we open the bed file and we shift the alignements
  my $bed_it = CracTools::Utils::bedFileIterator($bed_file); 
  while (my $bed_line = $bed_it->()) {
    my $new_line = { 
      chr     => $bed_line->{chr}, # TODO special treatment for fusions
      name    => $bed_line->{name},
      strand  => $bed_line->{strand},
      blocks  => [],
    };
    # We loop over each block of the bed alignement
    foreach my $block (@{$bed_line->{blocks}}) {
      # First we get shifted intervals
      my @shifted_intervals = @{$self->liftover->shiftInterval(
        $bed_line->{chr},
        $block->{ref_start},
        $block->{ref_end},
      )};
      foreach my $shifted_interval (@shifted_intervals) {
        if(!defined $new_line->{start}) {
          $new_line->{start}  = $shifted_interval->{start}; 
          $new_line->{chr}    = $shifted_interval->{chr};
          $new_line->{strand} = $shifted_interval->{strand};
        }
        my $block_start;
        # If this block is collinear with the bed_line genomic coordinate
        # we only update the block_start according to the line start
        if($new_line->{chr} eq $shifted_interval->{chr} &&
          $new_line->{strand} eq $shifted_interval->{strand}) {
          $block_start = $shifted_interval->{start} - $new_line->{start};
        # Otherwise it is a chimeric alignment and we have to put some
        # dirty things in the bed alignment
        } else {
          $block_start = $shifted_interval->{chr}."@".$shifted_interval->{strand}.$shifted_interval->{start};
        }
        # If the prev block is glued to this one (due to a deletion)
        # we merge them into one block.
        my $prev_block = @{$new_line->{blocks}}? $new_line->{blocks}->[$#{$new_line->{blocks}}]: undef;
        if(defined $prev_block && $prev_block->{start}+$prev_block->{size} == $block_start) {
          $prev_block->{size} += $shifted_interval->{end} - $shifted_interval->{start} + 1,
        # Otherwise we add a new block
        } else {
          push @{$new_line->{blocks}}, {
            size  => $shifted_interval->{end} - $shifted_interval->{start} + 1,
            start => $block_start,
          };
        }
      }
      # Then we get overlapping mutations
      map { push(@{$overlapped_mutations{$_}},$bed_line->{name}) } @{$self->mutation_query->fetchByRegion(
        $bed_line->{chr},
        $block->{ref_start},
        $block->{ref_end},
      )};
    }
    my $nb_blocks = @{$new_line->{blocks}};
    # Set end position of the alignement using the last block
    $new_line->{end} = $new_line->{start} + $new_line->{blocks}->[$nb_blocks-1]->{start} + $new_line->{blocks}->[$nb_blocks-1]->{size};
    CracTools::SimCT::Utils::printBEDLine($output_fh,$new_line);
    # TODO we should construct a SAM file instead
  }
  # Now we write a vcf file with expressed mutations
  if(defined $vcf_fh) {
    foreach my $mut (sort { $a->pos <=> $b->pos } keys %overlapped_mutations) {
      # Get raw VCF record
      my $vcf_line = $mut->getVCFRecord;
      # Update information with read names and depth
      $vcf_line->{id} = join(',',@{$overlapped_mutations{$mut}});
      $vcf_line->{info}->{DP} = [scalar @{$overlapped_mutations{$mut}}];
      CracTools::SimCT::Utils::printVCFLine($vcf_fh,$vcf_line);
    }
  }
}

1;

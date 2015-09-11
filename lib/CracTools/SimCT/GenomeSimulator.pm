package CracTools::SimCT::GenomeSimulator;
# ABSTRACT: Simulate a mutated genome by introducing indel, snv and fusion genes in a reference genome

use strict;
use warnings;

use List::Util qw(min max);
use Carp;
use Storable; # In order to clone the GenomeSimulator
use Data::Dumper;

use CracTools::Utils;
use CracTools::GenomeMask;
use CracTools::Interval::Query;
use CracTools::SimCT::Const;

=head2 new

=cut

sub new {
  my $class = shift;
  my %args = @_;

  # A hash with reference name as keys and fasta file as values
  my $reference_sequence_files = $args{reference_sequence_files};
  # GTF file with annotation over the genome
  my $annotation_file = $args{annotation_file};

  my $self = bless {
    reference_sequence_files => $reference_sequence_files,
    annotation_file => $annotation_file,
    genes                     => {},
    mutations                 => [], # Hold: insertions, deletions and substitutions
    fusions                   => [], # Hold: gene fusions
  }, $class;

  $self->_init();

  return $self;

}

sub _init {
  my $self = shift;

  # First we create a genomeMask over the genome
  # but we have to comptue the reference lengths first
  my %references_length;
  foreach my $reference ($self->references) {
    my $seq_it = CracTools::Utils::seqFileIterator($self->getReferenceFile($reference));
    my $entry = $seq_it->();
    $references_length{$reference} = length $entry->{seq};
  }
  $self->{references_length} = \%references_length;
  $self->{genome_mask} = CracTools::GenomeMask->new(genome => \%references_length);

  # Now we read the annotations and load them into memory
  my $gtf_it = CracTools::Utils::gffFileIterator($self->getAnnotationFile,'gtf');
  while(my $annot = $gtf_it->()) {
    # We only consider exon annotations
    next if $annot->{feature} ne 'exon';
    # We also do not load annotation that correspond to references we do not have
    next if !defined $self->getReferenceFile($annot->{chr});
    # Convert to 0-based
    $annot->{start}--;
    $annot->{end}--;
    my $gene_id = $annot->{attributes}->{gene_id};
    my $gene = $self->getGene($gene_id);
    # If this gene is already known, we adjust its boundaries
    # if needed
    if(defined $gene) {
      $gene->{start} = $annot->{start} if $annot->{start} < $gene->{start};
      $gene->{end} = $annot->{end} if $annot->{end} > $gene->{end};
    # otherwise we add a new gene to the dictionnary
    } else {
      $self->{genes}->{$gene_id} = {
        chr     => $annot->{chr},
        start   => $annot->{start},
        end     => $annot->{end},
        strand  => $annot->{strand},
        exons   => {},
      };
      $gene = $self->getGene($gene_id);
    }
    # Now we add the transcripts to its exons
    my $exon_key = _getExonKey($annot->{start},$annot->{end});
    push @{$gene->{exons}->{$exon_key}}, $annot->{attributes}->{transcript_id};
  }

  # We remove from the annotations, overlapping genes and exons
  my @gene_ids = sort { $self->getGene($a)->{chr} cmp $self->getGene($b)->{chr} || $self->getGene($a)->{start} <=> $self->getGene($b)->{start}} $self->genes;
  my $prev_gene_id = shift @gene_ids;
  foreach my $gene_id (@gene_ids) {
    my $prev_gene = $self->getGene($prev_gene_id);
    my $gene = $self->getGene($gene_id);
    if($gene->{chr} eq $prev_gene->{chr} && $gene->{start} <= $prev_gene->{end}) {
      $self->removeGene($gene_id);
    } else {
      # Sort exons by start position to remove overlapping ones
      my @exons_ids = sort { _getExonStart($a) <=> _getExonStart($b) } $self->{exons};
      my $prev_exon_id = shift @exons_ids;
      foreach my $exon_id (@exons_ids) {
        if(_getExonStart($exon_id) <= _getExonEnd($prev_exon_id)) {
          $self->removeExon($gene_id,$exon_id);
        } else {
          $prev_exon_id = $exon_id;
        }
      }
      $prev_gene_id = $gene_id;
    }
  }
}

# Return the annotation filename
sub getAnnotationFile {
  my $self = shift;
  return $self->{annotation_file};
}

# Return the genome mask
sub genomeMask {
  my $self = shift;
  return $self->{genome_mask};
}

# Retrun an array of all reference names
sub references {
  my $self = shift;
  return keys %{$self->{reference_sequence_files}};
}

# Given a reference name, return its filename
sub getReferenceFile {
  my $self = shift;
  my $reference = shift;
  return $self->{reference_sequence_files}->{$reference};
}
 
# Given a reference name, return its length
sub getReferenceLength {
  my $self = shift;
  my $reference = shift;
  return $self->{references_length}->{$reference};
}

# Return an array of all stored gene ids
sub genes {
  my $self = shift;
  return keys %{$self->{genes}};
}

# Given a gene_id, return the corresponding gene entry
sub getGene {
  my $self = shift;
  my $gene_id = shift;
  return $self->{genes}->{$gene_id};
}

# Given a gene_id, return its exon keys (see _getExonKey)
sub exons {
  my $self = shift;
  my $gene_id = shift;
  return keys %{$self->{genes}->{$gene_id}->{exons}} if defined $self->{genes}->{$gene_id};
}

# Given a gene_id, and an exon key, return all the transcripts_ids
sub transcripts {
  my $self = shift;
  my ($gene_id, $exon_key) = @_;
  return @{$self->getGene($gene_id)->{exons}->{$exon_key}} if defined $self->getGene($gene_id);
}

# Given a gene_id, remove it from the loaded annotations
sub removeGene {
  my $self = shift;
  my $gene_id = shift;
  delete $self->{genes}->{$gene_id} if defined $gene_id;
}

# Given a gene_id and a exon hash key, remove it from the loaded annotations
sub removeExon {
  my $self = shift;
  my ($gene_id, $exon_key) = @_;
  delete $self->getGene($gene_id)->{exons}->{$exon_key} if defined $gene_id && defined $self->getGene($gene_id) && defined $exon_key;
}

# Given a chr and pos and a nucleotid sequence, we mute the genome with an insertion
# Return true is the insertion was correctly added
sub addInsertion {
  my $self = shift;
  my ($chr, $pos, $inserted_sequence) = @_;
  if($pos >= $self->getReferenceLength($chr)) {
    carp "Insertion position is greater that the reference length";
  } elsif(!$self->genomeMask->getPos($chr,$pos)) {
    push @{$self->{mutations}}, {
      type  => 'ins',
      chr               => $chr, 
      pos               => $pos, 
      inserted_sequence => $inserted_sequence,
    };
    $self->genomeMask->setPos($chr,$pos);
    return 1;
  }
  return 0;
}

# Give a chr, a pos and a length, we mute the genome with a deletion
sub addDeletion {
  my $self = shift;
  my ($chr, $pos, $length) = @_;
  if($pos + $length > $self->getReferenceLength($chr)) {
    carp "Deletion longer than the reference";
  } elsif($self->genomeMask->getNbBitsSetInRegion($chr,$pos,$pos+$length-1) == 0) {
    push @{$self->{mutations}}, {
      type    => 'del',
      chr     => $chr,
      pos     => $pos,
      length  => $length,
    };
    $self->genomeMask->setRegion($chr,$pos,$pos+$length-1);
    return 1;
  }
  return 0;
}

# Given a chr, a pos and a nucleotide, we mute the genome with a substitution
sub addSubstitution {
  my $self = shift;
  my ($chr, $pos, $new_nuc) = @_;
  if($pos >= $self->getReferenceLength($chr)) {
    carp "Substitution position is greater that the reference length";
  } elsif(length $new_nuc != 1) {
    carp "Substitution is not composed of 1 nucleotid";
  } elsif(!$self->genomeMask->getPos($chr,$pos)) {
    push @{$self->{mutations}}, {
      type    => 'sub',
      chr     => $chr,
      pos     => $pos,
      new_nuc => $new_nuc,
    };
    $self->genomeMask->setPos($chr,$pos);
    return 1;
  }
  return 0;
}

# Give two gene_ids and two exon_ids
sub addFusion {
  my $self = shift;
  my ($gene_A_id, $exon_A_id, $gene_B_id, $exon_B_id) = @_;
  my $gene_A = $self->getGene($gene_A_id);
  my $gene_B = $self->getGene($gene_B_id);
  if(defined $gene_A && defined $gene_B) {
    push @{$self->{fusions}}, {
      gene_A  => $gene_A_id,
      exon_A  => $exon_A_id,
      gene_B  => $gene_B_id,
      exon_B  => $exon_B_id
    };
  }
}

sub fusions {
  my $self = shift;
  if(defined $self->{fusions}) {
    return @{$self->{fusions}};
  } else {
    return ();
  }
}

# Return an array with all the mutations
sub mutations {
  my $self = shift;
  return @{$self->{mutations}};
}


sub liftoverQuery {
  my $self = shift;
  $self->generateLiftover if !defined $self->{liftover_query};
  return $self->{liftover_query};
}

# Create a hash key from exon coordinates
sub _getExonKey {
  my ($start, $end) = @_;
  return "$start,$end";
}

# Return the end position of an exon hash key
sub _getExonEnd {
  my $exon_key = shift;
  return (split(",",$exon_key))[1];
}

# Return the start position of an exon hash key
sub _getExonStart {
  my $exon_key = shift;
  return (split(",",$exon_key))[0];
}

# Generate the simulated genome in its current state
# There is three possible output files:
# - The fasta geneome
# - the mutation file
# - the gtf file
# Return a copy of the GenomeSimulator that is a "frozen" copy of the simulated genome
# the original one, can still continue to evolve
sub generateGenome {
  my $self            = shift;
  my %args            = @_;
  my $fasta_dir       = $args{fasta_dir};
  my $mut_file        = $args{mutation_file};
  my $gtf_file        = $args{gtf_file};
  my $gtf_output_fh   = CracTools::Utils::getWritingFileHandle($gtf_file) if defined $gtf_file;
  my $mut_output_fh   = CracTools::Utils::getWritingFileHandle($mut_file) if defined $mut_file;
  
  # We generate the simulated genome for each chr
  foreach my $chr (sort $self->references) {
    my $fasta_output_fh = CracTools::Utils::getWritingFileHandle("$fasta_dir/$chr.fa") if defined $fasta_dir;
    my $fasta_input_it  = CracTools::Utils::seqFileIterator($self->getReferenceFile($chr));
    my $chr_entry       = $fasta_input_it->(); # This load the entire FASTA entry into memory
    my $remainder       = 0; # What's left in the current FASTA line
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome

    # We update the fusion genes that involves this chromosome and
    # store the corresponding sequences
    foreach my $fusion ($self->fusions) {
      foreach my $side (('A','B')) {
        my $fused_gene = $self->getGene($fusion->{"gene_$side"});
        my $fused_exon = $fusion->{"exon_$side"};

        # We only consider this gene if it is located on the current chromosom
        next if $fused_gene->{chr} ne $chr;

        # Depending on the side and the strand, we retrieve the appropriate sequence
        # and we place the annotations relative to this fusion into the annotation
        # structure
        my $fused_sequence; 
        my $fused_start;

        if(($fused_gene->{strand} eq '+' && $side eq 'A') || 
           ($fused_gene->{strand} eq '-' && $side eq 'B')) {
          my $fusion_length = _getExonEnd($fused_exon) - $fused_gene->{start} + 1;
          $fused_sequence   = substr($chr_entry->{seq},$fused_gene->{start},$fusion_length);
          $fused_start      = $fused_gene->{start};
        } else {
          my $fusion_length = $fused_gene->{end} - _getExonStart($fused_exon) + 1;
          $fused_sequence   = substr($chr_entry->{seq},_getExonStart($fused_exon),$fusion_length);
          $fused_start      = _getExonStart($fused_exon);
        }

        # We reverse the sequence if its originate from the reverse strand
        $fused_sequence  = CracTools::Utils::reverseComplement($fused_sequence) if $fused_gene->{strand} eq '-';

        $fusion->{"sequence_$side"} = $fused_sequence;
        $fusion->{"start_$side"}    = $fused_start;
      }
    }

    # Write the header of fasta output
    print $fasta_output_fh ">".$chr_entry->{name}."\n" if defined $fasta_output_fh;
    
    # We create an array with all exons that will be modified when applying all the mutations
    # to the reference genome
    # Get genes annotations sorted by position
    # This does not modify the original loaded annotations
    my @exons;
    my @genes_sorted = sort { $self->getGene($a)->{start} <=> $self->getGene($b)->{start} } grep { $self->getGene($_)->{chr} eq $chr} $self->genes;
    foreach my $gene_id (@genes_sorted) {
      my $gene          = $self->getGene($gene_id);
      my @exons_sorted  = sort { _getExonStart($a) <=> _getExonStart($b) } $self->exons($gene_id);
      foreach my $exon_key (@exons_sorted) {
        my @transcript_ids = $self->transcripts($gene_id,$exon_key);
        push @exons, {
          chr             => $gene->{chr},
          feature         => 'exon',
          start           => _getExonStart($exon_key),
          end             => _getExonEnd($exon_key),
          strand          => $gene->{strand},
          gene_id         => $gene_id,
          transcript_ids  => \@transcript_ids,
        };
      }
    }
    my $index_exon_start  = 0;
    my $index_exon_end    = 0;

    # Get sorted mutations for this chromosome
    my @mutations_sorted = sort {$a->{pos} <=> $b->{pos}} grep {$_->{chr} eq $chr} $self->mutations;
    foreach my $mut (@mutations_sorted) {
      # We read the sequence before the mutation and print it in
      # the output fasta
      my $frag    = substr $chr_entry->{seq}, 0, $mut->{pos}-$index, "";
      $index      = $mut->{pos};
      $remainder  = _printFASTA($fasta_output_fh,$frag,$remainder) if defined $fasta_output_fh;
      
      # We update the annotations according to the previous introduced mutations
      # Remember that ther is now overlap betwwen annotation, because we
      # have filtered them out in the _init() method.
      while($index_exon_start < @exons && $exons[$index_exon_start]->{start} < $mut->{pos}) {
        $exons[$index_exon_start]->{start} += $offset;  
        $index_exon_start++;
      }
      while($index_exon_end < @exons && $exons[$index_exon_end]->{end} < $mut->{pos}) {
        $exons[$index_exon_end]->{end} += $offset;  
        $index_exon_end++;
      }

      # MUTATIONS CASES LOOP
      # - We loop over the mutations cases and for each one of then
      #   we update the $index, the $offset and print sequence to the FASTA output if needed
      #
      # 1. INSERTION CASE
      # - If this mutation is an insertion we print the inserted sequence in place
      # - Offset is moving forward
      # - Index on the original sequence is not moved
      if($mut->{type} eq 'ins') {
        $remainder = _printFASTA($fasta_output_fh,$mut->{inserted_sequence},$remainder) if defined $fasta_output_fh;
        $offset += length $mut->{inserted_sequence};
      # 2. DELETION CASE
      # - If it is a deletion, we read the deleted sequence and remove it
      #   from the original reference
      # - Index is moving forward
      # - Offset is decremented by the deletion length
      } elsif ($mut->{type} eq 'del') {
        my $index_after_del = $index + $mut->{length};
        # If the deletion overlaps an annotation, we have to update its coordinates
        # Three cases are possible :
        # 1. Deletion start before the exon and stop in the middle
        # 2. Deletion start in the middle of the exon and stop after
        # 3. Deletion completely overlap the exon
        my $found_overlap = 1;
        while($found_overlap) { 
          $found_overlap = 0;
          if($index_exon_start < @exons && $index <= $exons[$index_exon_start]->{start}
            && $index_after_del > $exons[$index_exon_start]->{start}) {
            if($index_after_del > $exons[$index_exon_start]->{end}) {
              # The exon is completely overlaped by the deletion, we remove it
              $exons[$index_exon_start]->{deleted} = 1;
            } else {
              $exons[$index_exon_start]->{start}  = $index; 
              $exons[$index_exon_start]->{end}    = $index + ($exons[$index_exon_start]->{end} - $index_after_del);
            }
            $index_exon_start++;
            $found_overlap = 1;
          }
          if($index_exon_end < @exons && $index <= $exons[$index_exon_end]->{end}
            && $index_after_del > $exons[$index_exon_end]->{end}) {
            $exons[$index_exon_end]->{end}        = $index - 1;
            $index_exon_end++;
            $found_overlap = 1;
          }
        }
        $mut->{deleted_sequence} = substr $chr_entry->{seq},0,$mut->{length},"";
        $index  += $mut->{length};
        $offset -= $mut->{length};
      # 3. SUBSTITUTION CASE
      # - If this is a substitution, we read the old nuclotide and insert the new one instead
      # - Index is moving forward from 1 position
      } elsif ($mut->{type} eq 'sub') {
        my $old_nuc = substr $chr_entry->{seq},0,1,"";
        if($old_nuc eq $mut->{new_nuc}) {
          #carp "Substitution alternative nucleotid was identical to the reference";
        } else {
          $mut->{old_nuc} = $old_nuc;
        }
        $remainder = _printFASTA($fasta_output_fh,$mut->{new_nuc},$remainder) if defined $fasta_output_fh;
        $index++;
      } else {
        carp("Unknown mutation type ".$mut->{type});
      }
      #CracTools::Utils::writeSeq($filehandle,$format,$seq_name,$seq,$seq_qual)
    }

    # Now we print what's left of the reference
    _printFASTA($fasta_output_fh,$chr_entry->{seq},$remainder);

    # Switch the last exons
    while($index_exon_start < @exons) {
      $exons[$index_exon_start]->{start} += $offset;  
      $index_exon_start++;
    }
    while($index_exon_end < @exons) {
      $exons[$index_exon_end]->{end} += $offset;  
      $index_exon_end++;
    }

    # And print the exons to the GTF
    if(defined $gtf_output_fh) { _printGTF($gtf_output_fh,$_) foreach grep {!defined $_->{deleted}} @exons;}
  }

  # Now we print an extra FASTA file with fusions
  # TODO fused exons should be either merged into on new exon or separated by some intron,
  # otherwise, Flux Won't like it
  my $fasta_output_fh = CracTools::Utils::getWritingFileHandle("$fasta_dir/chrFusions.fa") if defined $fasta_dir;
  print $fasta_output_fh ">chrFusions\n" if defined $fasta_output_fh;
  my $remainder       = 0;
  my $fusion_id       = 0;
  my $chr_fusion_pos  = 0;
  foreach my $fusion ($self->fusions) {
    # First we verify that we have the sequences of both part of the fusions
    # otherwise, one of the chromosome involve in the fusion was not available
    next if !defined $fusion->{sequence_A} || ! defined $fusion->{sequence_B};
    # We print the fasta sequence corresponding to the fusion
    $remainder = _printFASTA($fasta_output_fh,$fusion->{sequence_A}.$fusion->{sequence_B},$remainder) if defined $fasta_output_fh;
    # Update the position of the fusion in the "Fusions" chr
    $fusion->{chr_fusion_pos} = $chr_fusion_pos;
    # No we update the annotations and print them to the GTF
    my $fusion_gene_id = "fusion_$fusion_id";
    my @exons;
    foreach my $side (('A','B')) {
      my $fused_gene_id = $fusion->{"gene_$side"};
      my $fused_gene    = $self->getGene($fused_gene_id);
      my $fused_exon    = $fusion->{"exon_$side"};
      my @sorted_exons;

      if($fused_gene->{strand} eq '+') {
        @sorted_exons = sort { _getExonStart($a) <=> _getExonStart($b) } $self->exons($fused_gene_id);
      } else {
        @sorted_exons = sort { _getExonStart($b) <=> _getExonStart($a) } $self->exons($fused_gene_id);
      }

      foreach my $exon (@sorted_exons) {
        if(($fused_gene->{strand} eq '+' && $side eq 'A') || 
           ($fused_gene->{strand} eq '-' && $side eq 'B')) {
          last if _getExonStart($exon) > _getExonStart($fused_exon);
        } else {
          last if _getExonStart($exon) < _getExonStart($fused_exon);
        }

        my $exon_start;
        my $exon_length = _getExonEnd($exon) - _getExonStart($exon) + 1; 
        if($fused_gene->{strand} eq '+') {
          $exon_start = $chr_fusion_pos + (_getExonStart($exon) - $fusion->{"start_$side"});
        } else {
          $exon_start = $chr_fusion_pos + $fusion->{"start_$side"} + length($fusion->{"sequence_$side"}) - 1 - _getExonEnd($exon);
        }

        my $exon_key    = _getExonKey($exon_start,$exon_start + $exon_length - 1);
        # If this is the first exon of side B, we merge it into a new exon that is made of the
        # fusion of the two fused exons
        if($exon eq $fused_exon && $side eq 'B') {
          my $new_exon = pop @exons;
          $new_exon->{end} += $exon_length;
          push @exons, $new_exon;
        } else {
          push @exons, {
            chr             => "Fusions",
            feature         => 'exon',
            start           => $exon_start,
            end             => $exon_start + $exon_length - 1,
            strand          => "+",
            gene_id         => $fusion_gene_id,
            transcript_ids  => ["$fusion_gene_id.1"],
          };
        }

      }
      # We update the chr fusion pos for the next exon
      $chr_fusion_pos += $exons[$#exons]->{end} + 1;
    }
    $fusion_id++;
    if(defined $gtf_output_fh) { _printGTF($gtf_output_fh,$_) foreach @exons;}
  }
  close($fasta_output_fh);

  # TODO return of copy of the genomeSimulator
  my $clone = $self->clone;
  $clone->generateLiftover;
  return $clone;
}

# Create a "frozen" clone of the genome simulator
sub clone {
  my $self = shift;
  my $copy;

  # Copy non-ref objects
  foreach my $key (keys %$self) {
    if(!ref $self->{$key}) {
      $copy->{$key} = $self->{$key};
    }
  }

  # Copy references
  $copy->{genes}      = Storable::dclone($self->{genes});
  $copy->{mutations}  = Storable::dclone($self->{mutations});
  $copy->{fusions}    = Storable::dclone($self->{fusions});
  $copy->{references_length}    = Storable::dclone($self->{references_length});
  $copy->{reference_sequence_files}    = Storable::dclone($self->{reference_sequence_files});
  $copy->{frozen}     = 1;

  bless $copy, ref $self;
  return $copy;
}

# Generate a liftover interval_query object based on the
# current mutations.
sub generateLiftover {
  my $self = shift;
  my $liftover_query  = CracTools::Interval::Query->new();

  # Create a liftover object
  foreach my $chr (sort $self->references) {
    my $prev_pos        = 0;
    my $index           = 0; # Index of the original FASTA sequence
    my $offset          = 0; # Offset difference between the original and the mutated genome

    my @mutations_sorted = sort {$a->{pos} <=> $b->{pos}} grep {$_->{chr} eq $chr} $self->mutations;
    foreach my $mut (@mutations_sorted) {
      next if $mut->{type} eq 'sub';
      $prev_pos   = $index;
      $index      = $mut->{pos} - $offset;
      $liftover_query->addInterval($chr,$prev_pos,$index-1,1,{chr => $chr, offset => $offset});
      if ($mut->{type} eq 'ins') {
        $prev_pos = $index;
        $index   += length $mut->{inserted_sequence};
        $offset  -= length $mut->{inserted_sequence};
        $liftover_query->addInterval($chr,$prev_pos,$index-1,1,undef);
      } elsif ($mut->{type} eq 'del') {
        $offset += $mut->{length};
      }
    }
    # Add the last interval
    $liftover_query->addInterval($chr,$index,$self->getReferenceLength($chr)-$offset,1,{chr => $chr, offset => $offset});
  }

  foreach my $fusion ($self->fusions) {
    my $start_fusion = $fusion->{chr_fusion_pos};
    foreach my $side (('A','B')) {
      my $fused_gene = $self->getGene($fusion->{"gene_$side"});
      my $offset_value;
      if($fused_gene->{strand} eq '+') {
        $offset_value = { 
          chr     => $fused_gene->{chr}, 
          offset  => $fused_gene->{start} - $start_fusion,
        };
      } else {
        $offset_value = { 
          chr             => $fused_gene->{chr}, 
          reverse_offset  => $fused_gene->{end} + $start_fusion,
        };
      }
      $liftover_query->addInterval(
        "Fusions",
        $start_fusion,
        $fusion->{chr_fusion_pos}+length($fusion->{"sequence_$side"}) - 1,
        1,
        $offset_value,
      );
      $start_fusion += length($fusion->{"sequence_$side"});
    }
  }

  $self->{liftover_query} = $liftover_query;
  return $liftover_query;
}

# Given an interval over the simulated genome, return
# the corresponding interval(s) over the reference genome
# The output is a reference array with entry of type :
# { chr => $chr, start => $start, end => $end, strand => $strand }
sub shiftInterval {
  my $self = shift;
  my ($chr,$start,$end,$strand) = @_;
  my @shifted_intervals;
  my ($shift_intervals,$shift_offsets) = $self->liftoverQuery->fetchByRegion($chr,$start,$end,$strand);
  # Now we soft the intervals
  my @sorted_shift_intervals = sort { $shift_intervals->[$a]->{start} <=> $shift_intervals->[$b]->{start} } (0..@{$shift_intervals}-1);
  foreach my $i (@sorted_shift_intervals) {
    if(defined $shift_offsets->[$i]) {
      my ($shifted_start, $shifted_end);
      my $offset      = $shift_offsets->[$i]->{offset};
      my $new_strand  = !defined $strand ? '+' : $strand;
      # If we have an offset, we simply shift coordinates
      if(defined $offset) {
        $shifted_start = $start + $offset;
        $shifted_end   = min($shift_intervals->[$i]->{end} + $offset, $end + $offset);
      } else {
        # We have a "reverse_offset", then all calculous are done
        # the opposite way
        $offset = $shift_offsets->[$i]->{reverse_offset};
        $shifted_end   = $offset - $start;
        $shifted_start = min($offset - $shift_intervals->[$i]->{start}, $offset - $end);
        $new_strand    = $new_strand eq '+' ? '-' : '+';
      }
      my $new_interval = {
        start   => $shifted_start,
        end     => $shifted_end,
        chr     => $shift_offsets->[$i]->{chr},
        strand  => $new_strand,
      };
      $start += $new_interval->{end} - $new_interval->{start} + 1;
      push @shifted_intervals, $new_interval;
    }
  }
  return \@shifted_intervals;
}


# Given bed that contains alignement,
# shit the coordinates to match the original genome
# If a VCF filename is given, also generate this output
sub shiftAlignements {
  my $self      = shift;
  my %args      = @_;

  my $bed_file          = $args{bed_file};
  my $output_file       = $args{output_file};
  my $vcf_file          = $args{vcf_file};
  my $output_fh         = CracTools::Utils::getWritingFileHandle($output_file);

  # Now we open the bed file and we shift the alignements
  my $bed_it = CracTools::Utils::bedFileIterator($bed_file); 
  while (my $bed_line = $bed_it->()) {
    my $new_line = { 
      chr     => $bed_line->{chr},
      name    => $bed_line->{name},
      strand  => $bed_line->{strand},
      blocks  => [],
    };
    #my $start_pos;
    foreach my $block (@{$bed_line->{blocks}}) {
      my @shifted_intervals = $self->shiftIntervals($bed_line->{chr},$block->{ref_start},$block->{ref_end});
      foreach my $shifted_interval (@shifted_intervals) {
        if(!defined $new_line->{start}) {
          $new_line->{start}  = $shifted_interval->{start}; 
          $new_line->{chr}    = $shifted_interval->{chr};
          $new_line->{strand} = $shifted_interval->{strand};
        }
        my $block_start;
        if($new_line->{chr} eq $shifted_interval->{chr} && $new_line->{strand} eq $shifted_interval->{strand}) {
          $block_start = $shifted_interval->{start} - $new_line->{start};
        } else {
          $block_start = $shifted_interval->{chr}."@".$shifted_interval->{strand}.$shifted_interval->{start};
        }
        push @{$new_line->{blocks}}, {
          size  => $shifted_interval->{end} - $shifted_interval->{start} + 1,
          start => $block_start,
        };
      }
    }
    my $nb_blocks = @{$new_line->{blocks}};
    # Set end position of the alignement using the last block
    $new_line->{end} = $new_line->{start} + $new_line->{blocks}->[$nb_blocks-1]->{start} + $new_line->{blocks}->[$nb_blocks-1]->{size};
    _printBED($output_fh,$new_line);
    # TODO we should construct a SAM file instead
  }
}

# PRINT METHODS:

# Write the mutated fasta files
# @ARG Handle of the output file
# @ARG String to be outtputted
# @ARG Number of characters already written on
#      the current line
# @RETURN The number of characters printed on the
#         last line
sub _printFASTA($$$) {
  my ($handle, $string, $remainder) = @_;
  my $current = 0;
  if ($remainder > 0) {
    $current = min($CracTools::SimCT::Const::FASTA_LINE_LENGTH-$remainder,length $string);
    print $handle substr($string, 0, $current);
    if ($current == length $string) {
      return length($string)+$remainder;
    } else {
      print $handle "\n";
    }
  }
  for (; $current+$CracTools::SimCT::Const::FASTA_LINE_LENGTH <= length $string; $current += $CracTools::SimCT::Const::FASTA_LINE_LENGTH) {
    print $handle substr($string, $current, $CracTools::SimCT::Const::FASTA_LINE_LENGTH),"\n";
  }
  if ($current < length($string)) {
    print $handle substr($string, $current, length($string)-$current);
    $remainder = length($string)-$current;
  } else {
    $remainder = 0;
  }
  return $remainder;
}

sub _printGTF($$) {
  my ($fh,$exon) = @_;
  foreach my $transcript_id (@{$exon->{transcript_ids}}) {
    print $fh join ("\t",
      $exon->{chr},       # seqname
      "GenomeSimulator",  # source
      "exon",             # feature
      $exon->{start} + 1, # start
      $exon->{end} + 1,   # end
      ".",                # score
      $exon->{strand},    # strand
      ".",                # frame
      join(" ",           # attributes
        'gene_id "'.$exon->{gene_id}.'";',
        'transcript_id "'.$transcript_id.'";',
      ),
    ), "\n";
  }
}

sub _printBED($$) {
  my ($fh,$bed) = @_;
  print $fh join("\t",
    $bed->{chr},
    $bed->{start},
    $bed->{end},
    $bed->{name},
    0,
    $bed->{strand},
    '.',
    '.',
    '0,0,0',
    @{$bed->{blocks}},
    join(",",map { $_->{size} } @{$bed->{blocks}}),
    join(",",map { $_->{start} } @{$bed->{blocks}}),
  ),"\n";
}

1;

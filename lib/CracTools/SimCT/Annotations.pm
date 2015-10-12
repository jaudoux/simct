package CracTools::SimCT::Annotations;
# ABSTRACT: A collection of annotations

use Moose;
use Carp;

use CracTools::SimCT::Utils;
use CracTools::SimCT::Genome;
use CracTools::SimCT::Annotations::Exon;
use CracTools::SimCT::Annotations::Gene;

=head1 SYNOPSIS

A collection of annotations (gene => [exon => transcript])

Annotations are 0-based.

=cut

# We store genes in a hash where the id is the uniq gene_id
has genes_hash => (
  traits => ['Hash'],
  is  => 'ro',
  isa => 'HashRef[CracTools::SimCT::Annotations::Gene]',
  default => sub { {} },
  handles => {
    allGenes    => 'values',
    allGeneIds  => 'keys',
    getGene     => 'get',
    nbGenes     => 'count',
    removeGene  => 'delete',
  },
);

sub sortedGenes {
  my $self = shift;
  return sort { $a->chr cmp $b->chr || $a->start <=> $b->start } $self->allGenes;
}

sub addGene($) {
  my $self = shift;
  my $gene = shift;
  $self->genes_hash->{$gene->id} = $gene;
}

=head2 loadGTF


  $annotations->loadGTF($gtf_file,$genome);

If a genome is provided, only annotations wich references is available
in the genome are loaded.

=cut

sub loadGTF {
  my ($self,$gtf_file,$genome) = @_;

  # Now we read the annotations and load them into memory
  my $gtf_it = CracTools::Utils::gffFileIterator($gtf_file,'gtf');
  while(my $annot = $gtf_it->()) {
    # We only consider exon annotations
    next if $annot->{feature} ne 'exon';

    # TODO implement this feature
    # We also do not load annotation that correspond to references we do not have
    next if defined $genome && !defined $genome->getReferenceFile($annot->{chr});

    my $gene_id = $annot->{attributes}->{gene_id};
    next if !defined $gene_id;

    my $gene = $self->getGene($gene_id);

    # If this gene is not defined yet, we create a new entry
    if(!defined $gene) {
      $gene = CracTools::SimCT::Annotations::Gene->new(
        chr     => $annot->{chr},
        strand  => $annot->{strand},
        id      => $gene_id,
      );
      $self->addGene($gene);
    }

    my $exon = CracTools::SimCT::Annotations::Exon->new(
      start     => $annot->{start} - 1,
      end       => $annot->{end} - 1,
      gene      => $gene,
      transcripts => [$annot->{attributes}->{transcript_id}],
    );
  }
}

sub saveGTF {
  my ($self, $gtf_file,$liftover) = @_;
  my $fh = CracTools::Utils::getWritingFileHandle($gtf_file);
  $self->appendGTF($fh,$liftover);
}

sub appendGTF {
  my ($self,$fh,$liftover) = @_;
  my @gtf_buffer = ();
  foreach my $gene ($self->sortedGenes) {
    if(@gtf_buffer > 0 && $gtf_buffer[$#gtf_buffer]->{start} < $gene->start) {
      map { CracTools::SimCT::Utils::printGTFLine($fh,$_) } sort { $a->{start} <=> $b->{start} } @gtf_buffer;
      @gtf_buffer = ();
    }
    #my %transcripts_gtf_lines = ();
    #my @sorted_transcript = ();
    foreach my $exon ($gene->sortedExons) {
      my ($chr,$start,$end,$strand) = ($exon->chr,$exon->start,$exon->end,$exon->strand);
      # If we have a liftover object, we shift the exon's coordinates
      if(defined $liftover) {
        my $shifted_exon = $liftover->shiftAnnotation($chr,$start,$end,$strand);
        # If this exon does not exists in the new
        # genome we move to the next exon
        next if !defined $shifted_exon;
        $chr    = $shifted_exon->{chr};
        $start  = $shifted_exon->{start};
        $end    = $shifted_exon->{end};
        $strand = $shifted_exon->{strand};
      }
      foreach my $transcript_id ($exon->allTranscripts) {
        my $gtf_line = {
          chr         => $chr,
          feature     => 'exon',
          start       => $start + 1,
          end         => $end + 1,
          strand      => $strand,
          attributes  => {
            'gene_id'       => $gene->id,
            'transcript_id' => $transcript_id,
          },
        };
        push @gtf_buffer, $gtf_line;
        #if(!defined $transcripts_gtf_lines{$transcript_id}) {
        #  push @sorted_transcript, $transcript_id;
        #}
        #push @{$transcripts_gtf_lines{$transcript_id}}, $gtf_line;
      }
    }
    # Print GTF line grouped by 'transrcipt id'
    #foreach my $transcript (@sorted_transcript) {
    #  foreach my $gtf_line (@{$transcripts_gtf_lines{$transcript}}) {
    #    CracTools::SimCT::Utils::printGTFLine($fh,$gtf_line);
    #  }
    #}
  }
  map { CracTools::SimCT::Utils::printGTFLine($fh,$_) } @gtf_buffer;
}

1;

__END__

=head1 METHODS

=head2 new

Create a new 'CracTools::SimCT::Annotations' object

=head2 addGene('CracTools::SimCT::Annotations::Gene')

Add a new gene to the annotations

=head2 getGene('gene_id') => 'CracTools::SimCT::Annotations::Gene'

Given a gene_id return the gene

=head2 removeGene('gene_id')

Given a gene_id remove this gene from the collection (if it exists)

=head2 nbGenes => 'Int'

Return the number of genes

=head2 allGenes => Array['CracTools::SimCT::Annotations::Gene']

Return all genes contained in the annotations

=head2 allGeneIds => Array['gene_id']

Return all gene ids contained in the annotations

=head2 sortedGenes => Array['CracTools::SimCT::Annotations::Gene']

Return all genes sorted by 'chr' and 'start' position

=head2 loadGTF('gtf_file')

Load annotations contained if the gtf_file passed in argument

=head2 saveGTF('gtf_file')

Save the annotations in the gtf_file passed in argument

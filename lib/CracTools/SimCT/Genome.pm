package CracTools::SimCT::Genome;
# ABSTRACT: A genome as a set of FASTA references (files)

use Moose;
use CracTools::Utils;

# Automatically set by reading the reference sequence files
has references_length => (
  traits => ['Hash'],
  is => 'ro',
  isa => 'HashRef[Int]',
  default => sub { {} },
  handles => {
    getReferenceLength => 'get',
  },
);

has reference_sequence_files => (
  traits => ['Hash'],
  is => 'ro',
  isa => 'HashRef[Str]',
  required => 1, 
  trigger => sub {
    my $self = shift;
    foreach my $reference ($self->references) {
      my $seq_it = CracTools::Utils::seqFileIterator($self->getReferenceFile($reference));
      my $entry = $seq_it->();
      $self->references_length->{$reference} = length $entry->{seq};
    }
  },
  handles => {
    getReferenceFile => 'get',
    references       => 'keys',
  },
);

sub sortedReferences {
  my $self = shift;
  return sort { $a cmp $b } $self->references;
}

sub getReferenceSeq($) {
  my $self = shift;
  my $ref  = shift;
  my $fasta_input_it  = CracTools::Utils::seqFileIterator($self->getReferenceFile($ref));
  return $fasta_input_it->()->{seq}; # This load the entire FASTA entry into memory
}

1;

__END__

=head1 METHODS

=head2 new

  Arg [reference_sequence_files] : 'HashRef[Files]' - A hash references
                                    where keys are reference names and
                                    values the reference FASTA files

Create a new 'CracTools::SimCT::Genome' object

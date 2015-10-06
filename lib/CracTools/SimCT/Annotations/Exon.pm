package CracTools::SimCT::Annotations::Exon;
# ABSTRACT: An Exon

use Moose;

use CracTools::SimCT::Annotations::Gene;

has gene  => (
  is        => 'ro',
  isa       => 'CracTools::SimCT::Annotations::Gene',
  required  => 1,
  trigger   => sub {
    my $self = shift;
    $self->gene->addExon($self);
  }
);

has start => (is => 'rw', isa => 'Int', required => 1);

has end   => (is => 'rw', 
  isa       => 'Int',
  required  => 1,
  # Verify that 'end' is greater than start
  trigger   => sub {
    my $self = shift;
    if($self->end < $self->start) {
      die "End position (".$self->end.") must be greater (or equal) to start (".$self->start.")";
    }
  },
);

has transcripts => (
  traits  => ['Array'],
  is      => 'ro',
  isa     => 'ArrayRef[Str]',
  default => sub { [] },
  handles => {
    addTranscript   => 'push',
    allTranscripts  => 'elements',
  },
);

# Use gene value for chr
sub chr {
  my $self = shift;
  return $self->gene->chr;
}

# Use gene value for strand
sub strand {
  my $self = shift;
  return $self->gene->strand;
}

sub length {
  my $self = shift;
  return $self->end - $self->start + 1;
}

# Remove a transcript given its id
sub removeTranscript($) {
  my $self = shift;
  my $transcript  = shift;
  my @transcripts = grep { $_ ne $transcript } @{$self->transcripts};
  $self->transcripts(\@transcripts); 
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 ACCESSORS

=head2 chr => 'String'

Getter for the exon chr (retrieved from the gene)

=head2 strand => ['+'|'-']

Getter for the exon strand (retrieved from the gene)

=head2 start => 'Int'

Getter/setter for the exon's start

=head2 end => 'Int'

Getter/setter for the exon's end. End should be greater or equal to start

=head1 METHODS

=head2 new

  Arg [start] : 'Int' - exon's start position
  Arg [end]   : 'Int' - exon's end position
  Arg [gene]  : 'CracTools::SimCT::Annotations::Gene' - exon's gene

Create a new 'CracTools::SimCT::Annotations::Exon' object

=head2 addTranscript('transcript_id')

Add a new transcript id for this exon

=head2 removeTranscript('transcript_id')

Remove a transcript_id given its name

=head2 allTranscripts() => Array['transcript_id']

Return all transcript ids

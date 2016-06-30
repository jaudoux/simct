package CracTools::SimCT::Utils;
# ABSTRACT: Utility subroutines and types

use Moose::Util::TypeConstraints;

use List::Util qw(min max);
use CracTools::SimCT::Const;

subtype 'Strand',
  as 'Str',
  where { $_ =~ /^(\+|-)$/ },
  message { "Strand must be either '+' or '-'" };

subtype 'DNA',
  as 'Str',
  where { $_ =~ /^[ACGNTRYSWKMBDHVatgcnryswkmbdhv]*$/ },
  message { "A DNA Sequence must be encoded over the alphabet defined by IUPAC code" };

subtype 'DNAnuc',
  as 'Str',
  where { $_ =~ /^[ACGNTRYSWKMBDHVatgcnryswkmbdhv]$/ },
  message { "A DNA nucleotide must be encoded over the alphabet defined by IUPAC code" };

subtype 'CigarOperator',
  as 'Str',
  where { $_ =~ /^[MIDNSHPX=]$/ },
  message { "A Cigar operator must respect the SAM format specifications"},

subtype 'Natural',
  as 'Int',
  where { $_ > 0 };

sub reverseStrand($) {
  my $strand = shift;
  $strand eq '+'? return '-' : return '+';
}

sub printGTFLine($$) {
  my ($fh,$gtf) = @_;
  my $line = join ("\t",
    $gtf->{chr},       # seqname
    defined $gtf->{source}? $gtf->{source} : "GenomeSimulator",  # source
    $gtf->{feature},   # feature
    $gtf->{start}, # start
    $gtf->{end},   # end
    defined $gtf->{score}? $gtf->{score} : ".",                # score
    $gtf->{strand},    # strand
    defined $gtf->{frame}? $gtf->{frame} : ".",                # frame
    join(" ",map { $_ .' "'.$gtf->{attributes}->{$_}.'";' } sort keys %{$gtf->{attributes}})
  )."\n";
  #print STDERR "$line";
  print $fh $line;
}

sub printVCFLine($$) {
  my ($fh,$vcf) = @_;
  print $fh join("\t",
    $vcf->{chr},
    $vcf->{pos},
    $vcf->{id},
    $vcf->{ref},
    $vcf->{alt},
    '.',
    'PASS',
    join(";",map { $_."=".join(',',@{$vcf->{info}->{$_}})  } sort keys %{$vcf->{info}}),
  ),"\n";
}

sub printBEDLine($$) {
  my ($fh,$bed) = @_;
  my $bed_line = join("\t",
    $bed->{chr},
    $bed->{start},
    $bed->{end},
    $bed->{name},
    defined $bed->{score}? $bed->{score} : 0,
    $bed->{strand}
  );
  # If bed line is in 12-field format we append blocks informations
  if(defined $bed->{blocks}) {
    $bed_line = join("\t",
      $bed_line,
      '.',
      '.',
      '0,0,0',
      join(",",map { $_->{size} } @{$bed->{blocks}}),
      join(",",map { $_->{start} } @{$bed->{blocks}}),
    );
  }
  print $fh $bed_line, "\n";
}

# Write the mutated fasta files
# @ARG Handle of the output file
# @ARG String to be outtputted
# @ARG Number of characters already written on
#      the current line
# @RETURN The number of characters printed on the
#         last line
sub printFASTA($$$) {
  my ($handle, $string, $remainder) = @_;
  my $line_length = $CracTools::SimCT::Const::FASTA_LINE_LENGTH;
  my $current = 0;
  # Print what's left in the line
  if ($remainder > 0) {
    $current = min($line_length-$remainder,length $string);
    print $handle substr($string, 0, $current);
    if ($current == length $string) {
      return length($string)+$remainder;
    } else {
      print $handle "\n";
    }
  }
  # Print one or more full line(s)
  for (; $current+$line_length <= length $string; $current += $line_length) {
    print $handle substr($string, $current, $line_length),"\n";
  }
  # Print what's left in a new line
  if ($current < length($string)) {
    print $handle substr($string, $current, length($string)-$current);
    $remainder = length($string)-$current;
  } else {
    $remainder = 0;
  }
  return $remainder;
}

1;

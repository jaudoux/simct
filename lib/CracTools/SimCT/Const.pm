use strict;
use warnings;
package CracTools::SimCT::Const;
# ABSTRACT: SimCT constants


# Tolerance on mapping position
our $SUB_RATE             = 0.005;
our $INS_RATE             = 0.001;
our $DEL_RATE             = 0.001;
our $MAX_INS              = 15;
our $MAX_DEL              = 15;
our $FASTA_LINE_LENGTH    = 60;
our $DEBUG                = 0;
our $FLUX_BINARY          = 'flux-simulator';
our $FLUX_OUTPUT_BASENAME = "fluxSimulator";
our $CHR_FUSIONS          = "Fusions";
our $MAX_SPLICE_LENGTH    = 300000;

1;

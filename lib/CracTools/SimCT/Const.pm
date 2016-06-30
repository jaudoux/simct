use strict;
use warnings;
package CracTools::SimCT::Const;
# ABSTRACT: SimCT constants

# Tolerance on mapping position
our $SUB_RATE             = 0.01;
our $INS_RATE             = 0.005;
our $DEL_RATE             = 0.005;
our $MAX_INS              = 15;
our $MAX_DEL              = 15;
our $FASTA_LINE_LENGTH    = 60;
our $DEBUG                = 0;
our $FLUX_BINARY          = 'flux-simulator';
our $FLUX_OUTPUT_BASENAME = "fluxSimulator";
our $CHR_FUSIONS          = "Fusions";
our $MAX_SPLICE_LENGTH    = 300000;
our $OUTPUT_DIRECTORY     = "simCT_simulation";

# Read simulator default parameters
our $UNIQ_IDS               = 1;
our $DISABLE_ERROR_ENCODING = 0;

# Flux-simulator default parameters
our $NB_MOLECULES         = 5000000;
our $NB_READS             = 1000000;
our $READS_LENGTH         = 100;
our $FRAGMENT_LENGTH      = 300;
our $FRAGMENT_SD          = 75;

1;

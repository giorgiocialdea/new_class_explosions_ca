#!/usr/bin/env python
# coding: utf-8
"""Global numerical parameters for the computer-assisted bounds.

Changing values here changes the precision, interval splitting, and amount of
printing from the code.
"""

from sage.all import RealBallField, RealIntervalField

############################
##### TYPES PRECISION ######
############################

# Number of bits used by Sage's real ball / interval arithmetic.
PRECISION_BITS = 200

# RBF is real ball arithmetic; RIF is real interval arithmetic.  The verifier
# uses these to keep rigorous lower/upper bounds instead of ordinary floats.
RBF = RealBallField(PRECISION_BITS)
RIF = RealIntervalField(PRECISION_BITS)

############################
###### VERIFICATION ########
############################

# The formulas are written in the variable b = sqrt(alpha).  The program proves
# the requested inequalities for b in this closed interval.
B_INTERVAL = ("0", "14/9")

# Initial cut points inside B_INTERVAL.  The proof starts by splitting the
# original interval in these subintervals.
CUTS = ("1/1024", "1/64", "1/8", "1/2")

# Safety limits for adaptive subdivision.  If a box cannot be proved before
# MAX_DEPTH, or if too many boxes are generated, the verification reports FAIL.
MAX_DEPTH = 60
MAX_BOXES = 1_000_000

# Print progress every N processed boxes.  Zero disables this unless VERBOSE
# later turns on the default VERBOSE_COUNTER_I value.
PROGRESS_EVERY = 0

############################
###### PRINTING INFO #######
############################

# Quantity of printed information:
#   0 prints only the lemma/final status,
#   1 also prints row-level progress and the bounds read from bounds.txt,
#   2 additionally prints every accepted subinterval and its local margins.
VERBOSE = 0

# Default progress interval used when verbose output is requested.
VERBOSE_COUNTER_I = 512

# Set True to print every accepted subinterval.
PRINT_SUBINTERVALS = False

# Display setting: number of digits to use if later printing
# format numerical values explicitly.
PRINT_DIGITS = 20

#!/usr/bin/env python
# coding: utf-8
"""Global numerical parameters for the computer-assisted bounds."""

from sage.all import RealBallField, RealIntervalField

############################
##### TYPES PRECISION ######
############################
PRECISION_BITS = 200
RBF = RealBallField(PRECISION_BITS)
RIF = RealIntervalField(PRECISION_BITS)

############################
###### VERIFICATION ########
############################
# The variable used by the program is b = sqrt(alpha).
B_INTERVAL = ("0", "14/9")
CUTS = ("1/1024", "1/64", "1/8", "1/2")

MAX_DEPTH = 60
MAX_BOXES = 1_000_000
PROGRESS_EVERY = 0

############################
###### PRINTING INFO #######
############################
# Quantity of printed information:
#   0 prints only the lemma/final status,
#   1 also prints row-level progress and the bounds read from bounds.txt,
#   2 additionally prints every accepted subinterval and its local margins.
VERBOSE = 0
VERBOSE_COUNTER_I = 512
PRINT_SUBINTERVALS = False
PRINT_DIGITS = 20

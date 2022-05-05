#!/bin/sh

# Make sure files from previous run are deleted
rm database.hmm
rm *.hmm.h3*

# Concat all files into a single database file
cat ./hmm/*.hmm >> database.hmm

# Press database file into a binary that hmmer can read
hmmpress database.hmm

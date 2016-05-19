#!/bin/bash

# --------------- HELP -------------------------
# First argument is the question to parameterize the solver with (optional if second argument ommitted, default is q3)
# Second argument is the input file (optional, default is data/2d_example.opp)
# File example for a 2D solver : 
# 4 
# 6
# 7
# 1 3 4
# 2 6 2
# 3 2 5
# 4 3 3
# where k = 4, n = 6, m = 7, X(1) = 3, Y(1) = 4, X(2) = 6, Y(2) = 2, etc.
# File example for a 3D solver (required for q6 and q7) :
# 4 
# 6
# 7
# 10
# 1 3 4 8
# 2 6 2 9
# 3 2 5 1
# 4 3 3 2
# where h = 10, Z(1) = 8, Z(2) = 9, etc.
# ----------------------------------------------

question="q3"
opp_file="data/2d_example.opp"

if [[ $# -eq 1 ]]; then
	question=$1
elif [[ $# -eq 2 ]]; then
	question=$1
	opp_file=$2
fi

cd solver/ && make && cat ../$opp_file | grep -v '^#' | ./solver $question && cd ..
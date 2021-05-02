#!/usr/bin/env bash

# Load Modules
module load iqtree

# Define variables
MSA=$1

iqtree -s $1 -m GTR+I+G -nt 10


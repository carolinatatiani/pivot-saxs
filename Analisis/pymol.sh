#!/bin/bash

prot1=`echo $1'.pdb'`
prot2=`echo $2'.pdb'`
#create the script
echo "load "$prot1 $prot2"
as cartoon
color red, ss h
color yellow, ss s
color green, ss l+ ''
align "$1", "$2"
" > script.pml

# GO!
pymol script.pml
rm script.pml

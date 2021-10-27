#!/bin/bash
pdb=$1
exp=$2

 ./saxs.e -pdb ${pdb}.pdb -exp ${exp}.dat

#rm slurm* fort.100  hbonds.dat exp.dat open.pdb pivot* *bib input-script.sh list.dat

echo 'Fim do job!!'

#!/bin/sh
awk '/vcondcoords1/ {print $2}' $1 > vcond.basis
awk '/vcondcoords2/ {print $2,$4,$5,$6}' $1 >> vcond.basis
awk '/vcondbasis/ {print}' $1| sed 's/vcondbasis//g' >> vcond.basis
awk '/vcondrnelnorb/ {print}' $1| sed 's/vcondrnelnorb//g' > vcond.input
echo 'pair-distribution function in MOs' >> vcond.input
awk '/vcondpdfmo/ {print}' $1| sed 's/vcondpdfmo//g' >> vcond.input
echo 'MOs in AOs' >> vcond.input
awk '/vcondvmoao/ {print}' $1| sed 's/vcondvmoao//g' >> vcond.input
#!/bin/sh
awk '/vcondcoords1/ {print $2}' $1 > basis
awk '/vcondcoords2/ {print $2,$4,$5,$6}' $1 >> basis
awk '/vcondbasis/ {print}' $1| sed 's/vcondbasis//g' >> basis
awk '/vcondrnelnorb/ {print}' $1| sed 's/vcondrnelnorb//g' > input
echo 'pair-distribution function in MOs' >> input
awk '/vcondpdfmo/ {print}' $1| sed 's/vcondpdfmo//g' >> input
echo 'MOs in AOs' >> input
awk '/vcondvmoao/ {print}' $1| sed 's/vcondvmoao//g' >> input
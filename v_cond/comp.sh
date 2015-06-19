#!/bin/sh
strdate=$(date)
sed -i s/"unknowndate"/"$strdate"/ vcndex.f90
ifort -c $1 getgss.f90
ifort -c $1 integ.f90
ifort -c $1 p1dens.f90
ifort -c $1 rdsamo.f90
ifort -c $1 valAO.f90
ifort -c $1 cond.f90
ifort -c $1 vbrain.f90
ifort -c $1 vcndex.f90
ifort -openmp -o vcond.exe getgss.o p1dens.o rdsamo.o integ.o valAO.o cond.o vbrain.o vcndex.o -L /sara/sw/mkl-11.0.2/composer_xe_2013.2.146/mkl/lib/intel64 -lmkl_rt -lmkl_lapack95_lp64 -L ../slatec/static/ -lslatec
sed -i s/"$strdate"/"unknowndate"/ vcndex.f90

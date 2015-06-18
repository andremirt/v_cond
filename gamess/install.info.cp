#!/bin/csh
#   compilation configuration for GAMESS
#   generated on mac161.few.vu.nl
#   generated at Wed Jan 26 15:31:44 CET 2011
setenv GMS_PATH            /home/andresch/research/vcond/gamess_exp
#         machine type
setenv GMS_TARGET          linux64
#         FORTRAN compiler setup
setenv GMS_FORTRAN         ifort
setenv GMS_IFORT_VERNO     12
#         mathematical library setup
setenv GMS_MATHLIB         mkl
setenv GMS_MATHLIB_PATH    /opt/intel/mkl/lib/intel64
setenv GMS_MKL_VERNO       10
#         parallel message passing model setup
setenv GMS_DDI_COMM        sockets

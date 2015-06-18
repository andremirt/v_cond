$ SET NOVERIFY
$ SET DEFAULT GAMESS:
$!
$! 18 Oct 1996 - MWS
$!
$!   Activate and compile a single GAMESS source code module.
$!   interactively - $ @COMP module_name
$!           batch - $ SUBMIT COMP/PARAM=module_name/NAME=module_name
$!
$!   Select arch from ALPHA, VAX
$!
$ arch = "ALPHA"
$!
$ module = P1
$ IF module.EQS."" THEN WRITE SYS$OUTPUT "No module name given"
$ IF module.EQS."" THEN EXIT
$!
$ IF F$MODE().EQS."BATCH" THEN SET PROCESS/NAME=CMP_'module'
$!
$!   The following modules contain VAX specific statements
$!   which must be activated before compilation.
$!
$ IF module.EQS."ALDECI" THEN GOTO ACT
$ IF module.EQS."CPHF"   THEN GOTO ACT
$ IF module.EQS."CPROHF" THEN GOTO ACT
$ IF module.EQS."GRD2A"  THEN GOTO ACT
$ IF module.EQS."GUGDGA" THEN GOTO ACT
$ IF module.EQS."GUGDGB" THEN GOTO ACT
$ IF module.EQS."GUGDM2" THEN GOTO ACT
$ IF module.EQS."GUGEM"  THEN GOTO ACT
$ IF module.EQS."GUGSRT" THEN GOTO ACT
$ IF module.EQS."GVB"    THEN GOTO ACT
$ IF module.EQS."HSS2A"  THEN GOTO ACT
$ IF module.EQS."INT2A"  THEN GOTO ACT
$ IF module.EQS."IOLIB"  THEN GOTO ACT
$ IF module.EQS."LAGRAN" THEN GOTO ACT
$ IF module.EQS."LOCAL"  THEN GOTO ACT
$ IF module.EQS."LOCPOL" THEN GOTO ACT
$ IF module.EQS."MCCAS"  THEN GOTO ACT
$ IF module.EQS."MCQDPT" THEN GOTO ACT
$ IF module.EQS."MCQUD"  THEN GOTO ACT
$ IF module.EQS."MCSCF"  THEN GOTO ACT
$ IF module.EQS."MCTWO"  THEN GOTO ACT
$ IF module.EQS."MOROKM" THEN GOTO ACT
$ IF module.EQS."MP2"    THEN GOTO ACT
$ IF module.EQS."MP2GRD" THEN GOTO ACT
$ IF module.EQS."MTHLIB" THEN GOTO ACT
$ IF module.EQS."ORDINT" THEN GOTO ACT
$ IF module.EQS."RHFUHF" THEN GOTO ACT
$ IF module.EQS."TDHF"   THEN GOTO ACT
$ IF module.EQS."TRANS"  THEN GOTO ACT
$ IF module.EQS."TRFDM2" THEN GOTO ACT
$ IF module.EQS."UNPORT" THEN GOTO ACT
$!
$!   Anything else is pure FORTRAN, just use it.
$!
$ source_fortran := [.SOURCE]'module'.SRC
$ GOTO COMP
$!
$!   Activate the VMS version - *VMS in IOLIB/UNPORT, *I32 elsewhere
$!
$ ACT:
$ ASSIGN/USER_MODE [.SOURCE]'module'.SRC SRCIN
$ ASSIGN/USER_MODE SCR:'module'.FOR      CODEOUT
$ ASSIGN/USER_MODE SYS$INPUT             ACTIN
$ ASSIGN/USER_MODE SYS$OUTPUT            ACTOUT
$ IF ((module.EQS."IOLIB") .OR. (module.EQS."UNPORT"))
$ THEN
$ RUN GAMESS:ACTVTE
*VMS
$ ELSE
$ RUN GAMESS:ACTVTE
*I32
$ ENDIF
$ source_fortran := SCR:'module'.FOR
$!
$!   Compile the source code
$!   GUGDRT seems to have trouble in the HPO compiler
$!
$ COMP:
$!
$!   Select the architecture dependent compiler options
$!
$!   scalar VAX tested using VAX FORTRAN 5.7-133 under VMS 5.5-1
$!
$ IF arch.EQS."VAX"
$ THEN
$    SET VERIFY
$    FORTRAN 'source_fortran' -
          /OBJECT=SCR:'module'.OBJ -
          /OPTIMIZE -
          /NOG_FLOATING -
          /CHECK=(NOBOUNDS,NOOVERFLOW,NOUNDERFLOW) -
          /DEBUG=(NOSYMBOLS,TRACEBACK) -
          /LIST=SCR:'module'.LIS -
          /CROSS_REFERENCE -
          /SHOW=(INCLUDE,MAP,NOPREPROCESSOR) -
          /STANDARD=(NOSEMANTIC,SOURCE_FORM,NOSYNTAX) -
          /WARNINGS=(NODECLARATIONS,NOGENERAL)
$    SET NOVERIFY
$ ENDIF
$!
$!   AXP tested using VAX FORTRAN 6.0 under openVMS 1.0
$!
$ IF arch.EQS."ALPHA"
$ THEN
$    SET VERIFY
$    FORTRAN 'source_fortran' -
          /OBJECT=SCR:'module'.OBJ -
          /OPTIMIZE -
          /FLOAT=IEEE_FLOAT -
          /ASSUME=(NOACCURACY_SENSITIVE,NODUMMY_ALIASES) -
          /CHECK=(NOBOUNDS,NOOVERFLOW,NOUNDERFLOW) -
          /DEBUG=(NOSYMBOLS,TRACEBACK) -
          /LIST=SCR:'module'.LIS -
          /CROSS_REFERENCE -
          /SHOW=(INCLUDE,MAP,NOPREPROCESSOR) -
          /STANDARD=(NOSEMANTIC,SOURCE_FORM,NOSYNTAX) -
          /WARNINGS=(ALIGNMENT,NODECLARATIONS,GENERAL,TRUNCATED_SOURCE)
$    SET NOVERIFY
$ ENDIF
$!
$!  Now put the object code in the library
$!
$ LIBRARY GAMESS:GAMESS.OLB SCR:'module'.OBJ
$ DELETE SCR:'module'.OBJ;*
$ IF F$SEARCH("SCR:''module'.FOR").NES."" THEN DELETE SCR:'module'.FOR;*
$ EXIT

$ SET NOVERIFY
$!
$! 24 May 1993 - MWS - link-edit GAMESS
$!   interactive - $ @LKED exefile_name
$!         batch - $ SUBMIT LKED/PARAM=exefile_name
$!
$! Choose architecture from ALPHA or VAX
$!
$ arch = "ALPHA"
$!
$ IF F$MODE().EQS."BATCH" THEN SET PROCESS/NAME=LKED_GAMESS
$ IF F$MODE().EQS."BATCH" THEN GOTO LINK
$!
$ IF P1.EQS."" THEN INQUIRE P1 "Enter name of executable image"
$ EXEFILE = P1 + ".EXE"
$ WRITE SYS$OUTPUT "The excutable file will be ",EXEFILE
$ INQUIRE OK "Is this OK"
$ OK = F$EDIT(OK,"TRIM, UPCASE")
$ IF OK.NES."OK" THEN EXIT
$!
$ LINK:
$ SET DEFAULT GAMESS:
$ IF P1.EQS."" THEN EXIT
$ SET RMS  /DISK  /BLOCK_COUNT=64  /BUFFER_COUNT=2
$!
$! Link a scalar version
$!
$ IF ((arch.EQS."ALPHA")  .OR.  (arch.EQS."VAX"))
$ THEN
$    SET VERIFY
$    LINK GAMESS.OLB/LIBRARY/INCLUDE=GAMESS -
          /EXECUTABLE='P1'.EXE /NOMAP/TRACEBACK
$    SET NOVERIFY
$ ENDIF
$ EXIT

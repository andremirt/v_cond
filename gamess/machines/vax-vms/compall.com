$! 28 Apr 1993 - MWS
$!
$!  This file will compile all of the GAMESS source code on a VAX.
$!
$!  The next two lines require customization for your site.
$!  The architecture should be chosen from ALPHA, VAX
$!        being sure to type capital letters.
$!  The queue name should be your short time limit, fast turnaround queue.
$!
$ arch = "ALPHA"
$ queue = "SHORT"
$!
$!   Create a new object code library
$!
$ SAY := WRITE SYS$OUTPUT
$ SET DEFAULT GAMESS:
$!
$ IF arch.EQS."ALPHA" 
$ THEN
$    SAY "Creating a new GAMESS.OLB object code library..."
$    LIBRARY/CREATE=(BLO:14000,GLO:800,MOD:800,HIS:20) GAMESS.OLB
$ ENDIF
$!
$ IF arch.EQS."VAX" 
$ THEN
$    SAY "Creating a new GAMESS.OLB object code library..."
$    LIBRARY/CREATE=(BLO:6000,GLO:500,MOD:500,HIS:20) GAMESS.OLB
$ ENDIF
$!
$!   Get a listing of all *.SRC files, and open it as DIRLIST.
$!
$ DIRECTORY/NOSIZE/NODATE  [.SOURCE]*.SRC; -
     /NOHEADING/NOTRAILING/OUTPUT=DIR.LIS
$ OPEN DIRLIST DIR.LIS/READ
$!
$!   Submit a batch job to compile each .SRC file
$!   (submit to short time limit queue, but at lower priority
$!   so as to avoid hogging the queue.)
$!
$ SAY "Submitting GAMESS (re)compilation to queue ''queue'..."
$ SRCLOOP:
$    READ/END_OF_FILE=ENDCMP DIRLIST module
$    I = F$LOCATE("]",module) + 1
$    J = F$LOCATE(".SRC;",module) - I
$    module = F$EXTRACT(I,J,module)
$    IF (module.EQS."VECTOR") THEN GOTO SRCLOOP
$    SUBMIT COMP.COM /PRIORITY=1 /NOPRINTER /NONOTIFY /QUEUE='queue' -
            /PARAM='module' /NAME='module' /LOG='module'.LOG
$ GOTO SRCLOOP
$!
$ ENDCMP:
$ CLOSE DIRLIST
$ DELETE DIR.LIS;*
$ EXIT

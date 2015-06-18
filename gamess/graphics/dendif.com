$ SET NOVERIFY
$ IF F$MODE().NES."INTERACTIVE" THEN EXIT
$ SAY := WRITE SYS$OUTPUT
$!
$! 26 Apr 1993 - MWS
$!
$!   Command file to execute DENDIF program
$!
$!   This .COM file must assign the logical names:
$!      DDFINP  is to be the user's input file.
$!      DDFLOG  is the print output.
$!      FOR005  is the keyboard
$!      FOR006  is the screen
$!      FOR020  is the the molecular grid file for total density.
$!      FOR021,... are any grids to be subtracted from the total.
$!
$!
$ deck = P1
$ IF deck.EQS.""
$ THEN
$    SAY ""
$    SAY "Your input for DENDIF must be in a XXX.DDF file."
$    SAY "Your total density grid file from PLTORB must be SCR:XXX.RHO"
$    SAY ""
$    INQUIRE deck "What is XXX (EXIT to quit)"
$    IF deck.EQS."EXIT" THEN EXIT
$ ENDIF
$!
$ SAY ""
$ SAY "You will now be asked for the name(s) of the SCR:YYY.RHO grid"
$ SAY "files to be subtracted from the total density in SCR:''deck'.RHO"
$ SAY "Enter DONE when all (if any) subtraction grid sets have been"
$ SAY "given, EXIT to abort."
$ SAY ""
$ minus = 20
$LOOP:
$ minus=minus+1
$ INQUIRE subtract "What is the name of the next file's YYY (or DONE)"
$ IF subtract.EQS."EXIT" THEN EXIT
$ IF subtract.EQS."DONE" THEN GOTO GOTDDF
$ ASSIGN SCR:'subtract'.RHO FOR0'minus'
$ GOTO LOOP
$!
$GOTDDF:
$ minus=minus-21                  
$ SAY "''minus' sets of MO grid files will be subtracted."
$!
$ output = P2
$ IF output.EQS."" 
$ THEN
$    INQUIRE output "Do you want X-windows or PostScript output? (XW or PS)"
$ ENDIF
$!
$ SAY ""
$ SAY "Your printout will be in ''deck'.LIS."
$ SAY ""
$!
$ ASSIGN     'deck'.DDF DDFINP
$ ASSIGN     'deck'.LIS DDFLOG
$ ASSIGN            TT: FOR005
$ ASSIGN            TT: FOR006
$ ASSIGN SCR:'deck'.RHO FOR020
$!
$!   Run using DECwindows
$!
$ IF output.EQS."XW"
$ THEN
$    DEFINE MACX FALSE
$    RUN PLT:DENDIF_XW.EXE
$    DEASSIGN MACX
$ ENDIF
$!
$!   Run using PostScript
$!
$ IF output.EQS."PS"
$ THEN
$    ASSIGN 'deck'.PS  PSTNAM
$    RUN PLT:DENDIF_PS.EXE
$    SAY "Your PostScript image is now in file ''deck'.PS"
$    DEASSIGN PSTNAM
$ ENDIF
$!
$ DEASSIGN DDFINP
$ DEASSIGN DDFLOG
$ DEASSIGN FOR005
$ DEASSIGN FOR006
$ DEASSIGN FOR020
$ minus=minus+20
$ mm=20
$DEASS:
$ mm=mm+1
$ IF mm.gt.minus THEN EXIT
$ DEASSIGN FOR0'mm'
$ GOTO DEASS
$ EXIT

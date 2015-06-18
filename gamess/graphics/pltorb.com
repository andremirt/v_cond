$ SET NOVERIFY
$ SAY := WRITE SYS$OUTPUT
$!
$! 26 Apr 1993 - MWS
$!
$!   Command file to execute PLTORB program
$!
$!   Run interactively by $ @PLTORB, optionally giving three arguements
$!   which are the input file name, XW or PS, and Y or N to save the
$!   .RHO file for DENDIF.  Omitted args will be prompted for.
$!
$!   Run in batch by $ SUBMIT PLTORB/PARAM=(filename,PS only,Y or N),
$!   being sure to supply all three arguments.  This part could use
$!   some work, as it assumes inputs are in the home directory!
$!
$!   This .COM file must assign the logical names:
$!      PLTORB  is to be the user's input file.
$!      PLTLOG  is the print output
$!      PLTVEC  is the input formatted molecular orbitals.
$!      PRGRID  is the primitive AO's grid file.
$!      COGRID  is the contracted AO's grid file.
$!      MOGRID  is the molecular orbital grid file.
$!
$ deck = P1
$ IF deck.EQS.""
$ THEN
$    SAY "Your MO vectors (if any) should be in XXX.VEC, and"
$    SAY "your PLTORB card input should be in XXX.ORB."
$    SAY ""
$    INQUIRE deck "What is XXX?  (enter EXIT to quit)"
$    IF deck.EQS."EXIT" THEN EXIT
$ ENDIF
$!
$ output = P2
$ IF output.EQS."" 
$ THEN
$    INQUIRE output "Do you want X-windows or PostScript output? (XW or PS)"
$ ENDIF
$!
$ saveden = P3
$ IF saveden.EQS.""
$ THEN
$    SAY "It is possible to save the MO grid data for later use by DENDIF."
$    INQUIRE saveden "Should I save this file for DENDIF (Y/N)"
$ ENDIF
$!
$ IF saveden.EQS."Y" THEN SAY "Saving MO grids in file SCR:''deck'.RHO"
$ IF saveden.NES."Y" THEN SAY "MO grids will be discarded when done."
$!
$ SAY ""
$ SAY "Your printout will be in ''deck'.LIS."
$ SAY ""
$ ASSIGN            TT: FOR005
$ ASSIGN            TT: FOR006
$ ASSIGN     'deck'.ORB PLTORB
$ ASSIGN     'deck'.LIS PLTLOG
$ ASSIGN     'deck'.VEC PLTVEC
$ ASSIGN SCR:'deck'.F08 PRGRID
$ ASSIGN SCR:'deck'.F09 COGRID
$ IF saveden.NES."Y" THEN ASSIGN SCR:'deck'.F10 MOGRID
$ IF saveden.EQS."Y" THEN ASSIGN SCR:'deck'.RHO MOGRID
$!
$!   Run using DECwindows
$!
$ IF output.EQS."XW"
$ THEN
$    IF F$MODE().NES."INTERACTIVE" THEN EXIT
$    DEFINE MACX FALSE
$    RUN PLT:PLTORB_XW.EXE
$    DEASSIGN MACX
$ ENDIF
$!
$!   Run using PostScript
$!
$ IF output.EQS."PS"
$ THEN
$    ASSIGN 'deck'.PS PSTNAM
$    RUN PLT:PLTORB_PS.EXE
$    SAY "Your PostScript image is in the file ''deck'.PS"
$    DEASSIGN PSTNAM
$ ENDIF
$!
$ DEASSIGN FOR005
$ DEASSIGN FOR006
$ DEASSIGN PLTORB
$ DEASSIGN PLTLOG
$ DEASSIGN PLTVEC
$ DEASSIGN PRGRID
$ DEASSIGN COGRID
$ DEASSIGN MOGRID
$ DELETE SCR:'deck'.F*;*

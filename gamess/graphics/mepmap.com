$ SET NOVERIFY
$ IF F$MODE().NES."INTERACTIVE" THEN EXIT
$ SAY := WRITE SYS$OUTPUT
$!
$! 26 Apr 1993 - MWS
$!
$!   Command file to execute MEPMAP program
$!
$!   This .COM file must assign the logical names:
$!      MEPMAP is to be the MEPMAP input deck, in some disk file.
$!      MEPLOG is the print output.
$!      FOR005 is to be the user's terminal, for interactive input.
$!      FOR006 is the users's terminal, for interactive output.
$!
$ deck = P1
$ IF deck.EQS.""
$ THEN
$    INQUIRE deck "What is XXX of your XXX.MEP input file (or EXIT)"
$    IF deck.EQS."EXIT" THEN EXIT
$ ENDIF
$!
$ output = P2
$ IF output.EQS."" 
$ THEN
$    INQUIRE output "Do you want X-windows or PostScript output? (XW or PS)"
$ ENDIF
$!
$ ASSIGN 'deck'.MEP MEPMAP
$ ASSIGN 'deck'.LIS MEPLOG
$ ASSIGN        TT: FOR005
$ ASSIGN        TT: FOR006
$!
$!   Run using DECwindows
$!
$ IF output.EQS."XW"
$ THEN
$    DEFINE MACX FALSE
$    RUN PLT:MEPMAP_XW.EXE
$    DEASSIGN MACX
$ ENDIF
$!
$!   Run using PostScript
$!
$ IF output.EQS."PS"
$ THEN
$    ASSIGN 'deck'.PS PSTNAM
$    RUN PLT:MEPMAP_PS.EXE
$    SAY "Your PostScript image is now in file ''deck'.PS"
$    DEASSIGN PSTNAM
$ ENDIF
$!
$ DEASSIGN MEPMAP
$ DEASSIGN MEPLOG
$ DEASSIGN FOR005
$ DEASSIGN FOR006

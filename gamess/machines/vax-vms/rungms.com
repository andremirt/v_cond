$! 18 Mar 1997 - MWS - Run GAMESS under VMS.
$!
$!    Submit this job by $ SUBMIT RUNGMS/PARAM=xxx/NAME=xxx
$!
$ deck    = P1     ! name of .INP input deck.
$ exefile = P2
$ IF    deck.EQS."" THEN EXIT
$ IF exefile.EQS."" THEN exefile = "GAMESS"
$!
$!  Pick up a copy of the input.
$!  Remove the directory specifier to use this for production
$!
$ COPY GAMESS:[.TESTS]'deck'.INP SCR:'deck'.F05
$!
$ ASSIGN SCR:'deck'.IRC  IRCDATA
$ ASSIGN SCR:'deck'.F05  INPUT
$ ASSIGN SYS$OUTPUT      OUTPUT
$ ASSIGN SCR:'deck'.DAT  PUNCH
$ ASSIGN SCR:'deck'.F08  AOINTS
$ ASSIGN SCR:'deck'.F09  MOINTS
$ ASSIGN SCR:'deck'.F10  DICTNRY
$ ASSIGN SCR:'deck'.F11  DRTFILE
$ ASSIGN SCR:'deck'.F12  CIVECTR
$ ASSIGN SCR:'deck'.F13  NTNFMLA
$ ASSIGN SCR:'deck'.F14  CIINTS
$ ASSIGN SCR:'deck'.F15  WORK15
$ ASSIGN SCR:'deck'.F16  WORK16
$ ASSIGN SCR:'deck'.F17  CSFSAVE
$ ASSIGN SCR:'deck'.F18  FOCKDER
$ ASSIGN SCR:'deck'.F20  DASORT
$ ASSIGN SCR:'deck'.F23  JKFILE
$ ASSIGN SCR:'deck'.F24  ORDINT
$ ASSIGN SCR:'deck'.F25  EFPIND
$ ASSIGN SCR:'deck'.F26  PCMDATA
$ ASSIGN SCR:'deck'.F27  PCMINTS
$ ASSIGN SCR:'deck'.F30  DAFL30
$ ASSIGN SCR:'deck'.F50  MCQD50
$ ASSIGN SCR:'deck'.F51  MCQD51
$ ASSIGN SCR:'deck'.F52  MCQD52
$ ASSIGN SCR:'deck'.F53  MCQD53
$ ASSIGN SCR:'deck'.F54  MCQD54
$ ASSIGN SCR:'deck'.F55  MCQD55
$ ASSIGN SCR:'deck'.F56  MCQD56
$ ASSIGN SCR:'deck'.F57  MCQD57
$ ASSIGN SCR:'deck'.F58  MCQD58
$ ASSIGN SCR:'deck'.F59  MCQD59
$ ASSIGN SCR:'deck'.F60  MCQD60
$ ASSIGN SCR:'deck'.F61  MCQD61
$ ASSIGN SCR:'deck'.F62  MCQD62
$ ASSIGN SCR:'deck'.F63  MCQD63
$ ASSIGN SCR:'deck'.F64  MCQD64
$!
$ SET RMS_DEFAULT/DISK /BLOCK_COUNT=64 /BUFFER_COUNT=1
$ SET PROCESS/NAME=GMS_'deck'
$ RUN GAMESS:'exefile'.EXE
$!
$ DIRECTORY/SIZE=ALL/DATE SCR:'deck'.*
$ DELETE SCR:'deck'.F*;*
$!
$! NOTE - These last two destroy data from the test runs that
$!        you'd probably want to keep during production runs.
$!
$ DELETE SCR:'deck'.DAT;*
$ IF F$SEARCH("SCR:''deck'.IRC").NES."" THEN DELETE SCR:'deck'.IRC;*
$ EXIT

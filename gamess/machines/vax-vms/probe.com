$ SET NOVERIFY
$ SET DEFAULT GAMESS:
$ SAY := WRITE SYS$OUTPUT
$!
$! 01 Apr 1990 - MWS
$!
$!   This command file PROBEs the GAMESS program
$!   Written April, 1985 by Mike Schmidt
$!
$ SAY "PROBE options are:"
$ SAY "0. Quit"
$ SAY "1. List all changes to GAMESS source (history)"
$ SAY "2. List machine dependent source statements"
$ SAY "3. Alphabetical subprogram list"
$ SAY "4. Module list, showing all subprograms"
$ SAY "5. Subroutine cross reference from LIBRARIAN"
$ SAY "6. List most recent change to each module."
$ SAY "7. List most recent change, and count source lines."
$!-- SAY "8. All COMMON block occurences (sorted)"
$ INQUIRE OPT "Enter number of desired choice"
$ SAY ""
$ IF OPT.EQ.0 THEN EXIT
$ IF OPT.EQ.1 THEN GOTO HIST
$ IF OPT.EQ.2 THEN GOTO SRCH
$ IF OPT.EQ.3 THEN GOTO ROUT
$ IF OPT.EQ.4 THEN GOTO MODU
$ IF OPT.EQ.5 THEN GOTO XREF
$ IF OPT.EQ.6 THEN GOTO VERS
$ IF OPT.EQ.7 THEN GOTO VERS
$ SAY "unimplemented option"
$ EXIT
$! ------------------------------ history -------------------
$ HIST:
$!
$! Feb 26, 1985 - MWS - prepare a history of GAMESS source code
$!                      modifications from the comments at the
$!                      start of each source module.
$!
$!   Open the output HISTORY file.
$!
$ OPEN HISTORY SCR:HISTORY.LIS/WRITE
$ WRITE HISTORY "GAMESS HISTORY LISTING"
$ WRITE HISTORY "------ ------- -------"
$!
$!   Get a listing of all .SRC files, and open it as DIRLIST.
$!
$ DIRECTORY [.SOURCE]*.SRC /COLUMNS=1/OUTPUT=SCR:DIR.LIS
$ OPEN DIRLIST SCR:DIR.LIS/READ
$ READ DIRLIST FILESPEC
$ READ DIRLIST FILESPEC
$ READ DIRLIST FILESPEC  ! now positioned to first filename
$!
$!   Loop over each .SRC file
$!
$FILELOOP:
$ READ DIRLIST FILESPEC
$ IF FILESPEC.EQS."" THEN GOTO END
$ WRITE HISTORY " "
$ WRITE HISTORY " "
$ WRITE HISTORY "HISTORY OF MODULE ",FILESPEC
$ WRITE HISTORY " "
$ SAY "Reading module ",FILESPEC
$ OPEN SRCFILE 'FILESPEC'/READ
$!
$!   Read the comments at the top of this source file,
$!   stopping at the C*MODULE card, and write them to HISTORY.
$!
$CODELOOP:
$ READ SRCFILE CODELINE
$ FIRST2 = F$EXTRACT(0,2,CODELINE)
$ IF FIRST2.NES."C " THEN CLOSE SRCFILE
$ IF FIRST2.NES."C " THEN GOTO FILELOOP
$ WRITE HISTORY CODELINE
$ GOTO CODELOOP
$!
$!   HISTORY file is now complete
$!
$END:
$ SAY "GAMESS history is now in SCR:HISTORY.LIS"
$ CLOSE HISTORY
$ CLOSE DIRLIST
$ DELETE SCR:DIR.LIS;*
$ EXIT
$! ------------------------------ search --------------------
$ SRCH:
$!
$! Apr 13, 1985  - MWS - Search out machine dependent code lines
$!                       in all GAMESS source code module.
$!
$ SAY "Machine dependent code lines will be found,"
$ SAY "and written into a file for later inspection."
$ SAY "Current machine types are"
$ SAY "IBM, VAX, FPS, CRY, UNX, ETA, or ALL of the above."
$ INQUIRE MACH "Enter name of desired machine"
$ SAY "Be patient, the search takes a while..."
$ IF MACH.EQS."IBM" THEN SEARCH [.SOURCE]*.SRC "*IBM"/OUTPUT=SCR:IBM.LIS
$ IF MACH.EQS."VAX" THEN SEARCH [.SOURCE]*.SRC "*VAX"/OUTPUT=SCR:VAX.LIS
$ IF MACH.EQS."FPS" THEN SEARCH [.SOURCE]*.SRC "*FPS"/OUTPUT=SCR:FPS.LIS
$ IF MACH.EQS."CRY" THEN SEARCH [.SOURCE]*.SRC "*CRY"/OUTPUT=SCR:CRY.LIS
$ IF MACH.EQS."UNX" THEN SEARCH [.SOURCE]*.SRC "*UNX"/OUTPUT=SCR:UNX.LIS
$ IF MACH.EQS."ETA" THEN SEARCH [.SOURCE]*.SRC "*ETA"/OUTPUT=SCR:ETA.LIS
$ IF MACH.EQS."ALL" THEN SEARCH [.SOURCE]*.SRC /MATCH=OR -
          "*VAX","*IBM","*FPS","*CRY","*UNX","*ETA" /OUTPUT=SCR:ALL.LIS
$ SAY MACH," specific code is listed in SCR:",MACH,".LIS"
$ EXIT
$! ------------------------------ alphabetical routine list -----
$ ROUT:
$!
$! Apr 15, 1985 - MWS - list all GAMESS subprograms in alphabetical
$!                      order, and give module where they occur
$!
$ SAY "Preparing subprogram list, be patient..."
$ SEARCH [.SOURCE]*.SRC "C*MODULE" /EXACT/OUTPUT=MODULE.TMP -
         /NOHEADING/NOLOG/NUMBERS
$ SORT MODULE.TMP/KEY=(POSITION:31,SIZE:6,CHARACTER,ASCENDING) -
       SCR:ROUTINE.LIS/PROCESS=RECORD
$ DELETE MODULE.TMP;*
$ SAY "Alphabetical subprogram list is in file SCR:ROUTINE.LIS"
$ SAY "This file gives modules, and line numbers,"
$ SAY "where each subprogram can be found."
$ EXIT
$! ------------------------------ list modules, and subprograms ------
$ MODU:
$!
$! Apr 15, 1985 - MWS - list all GAMESS modules in alphabetical
$!                      order, and subprograms they contain.
$!
$ SAY "Preparing module list, be patient..."
$ SEARCH [.SOURCE]*.SRC "C*MODULE" /EXACT/OUTPUT=SCR:MODULE.LIS -
                          /HEADING/LOG/NUMBERS
$ SAY "Alphabetical module list is in file SCR:MODULE.LIS"
$ SAY "This file gives modules, the subprograms in each,"
$ SAY "and line numbers where each suprogram starts."
$ EXIT
$! -------------------------- Subroutine cross reference ------
$!
$! Apr 15, 1985 - MWS - Provide cross reference listing for GAMESS
$!
$ XREF:
$ SAY "Cross reference options"
$ SAY "0. Quit"
$ SAY "1. Backwards, routine X is called by..."
$ SAY "2. Forwards, routine X calls..."
$ INQUIRE OPT "Enter number of desired choice"
$ IF OPT.EQ.1 THEN GOTO BACK
$ IF OPT.EQ.2 THEN GOTO FORW
$ EXIT
$ BACK:
$ LIBRARY GAMESS.OLB/CROSS_REFERENCE=SYMBOL/OUTPUT=SCR:XREFBACK.LIS
$ SAY "Backwards cross reference is in SCR:XREFBACK.LIS"
$ EXIT
$ FORW:
$ LIBRARY GAMESS.OLB/CROSS_REFERENCE=MODULE/OUTPUT=SCR:XREFFOR.LIS
$ SAY "Forwards cross reference is in SCR:XREFFOR.LIS"
$ EXIT
$! ------------------------------ version -------------------
$ VERS:
$ SOURCE_LINES = 0
$ MODULES = 0
$ SAY "What you will see on the screen is also going into a disk file."
$ OPEN VERSION SCR:VERSION.LIS/WRITE
$ WRITE VERSION "module lines last changed by   comments"
$!
$!   Get a listing of all .SRC files, and open it as DIRLIST.
$!
$ DIRECTORY [.SOURCE]*.SRC; /COLUMNS=1/NOHEAD/NOTRAIL/OUTPUT=SCR:DIR.LIS
$ OPEN DIRLIST SCR:DIR.LIS/READ
$!
$!   Loop over each .SRC file
$!
$VERS_FILELOOP:
$ READ/END=VERS_ALLDONE DIRLIST FILESPEC
$ OPEN SRCFILE 'FILESPEC'/READ
$ MODULES = MODULES + 1
$ TEXT_LINES = 0
$!
$!  extract the name of the module from the full file name.
$!
$ LOCMOD = F$LOCATE("]",FILESPEC) + 1
$ LENMOD = F$LOCATE(".SRC;",FILESPEC) - LOCMOD
$ MOD = F$EXTRACT(LOCMOD,LENMOD,FILESPEC) + "      "
$ MODULE = F$EXTRACT(0,6,MOD)
$!
$!  grab the comment in the very first line
$!
$ READ SRCFILE VERSION_INFO
$ COMMENT = F$EXTRACT(1,F$LENGTH(VERSION_INFO)-1,VERSION_INFO)
$ IF OPT.EQ.6 THEN GOTO VERS_EOF_MODULE
$ TEXT_LINES = TEXT_LINES + 1
$!
$!  count the lines of code.  this is slow!
$!
$VERS_COUNT:
$ READ/END=VERS_EOF_MODULE SRCFILE TEXT
$ TEXT_LINES = TEXT_LINES + 1
$ GOTO VERS_COUNT
$!
$VERS_EOF_MODULE:
$ CLOSE SRCFILE
$ VERSION_TEXT = F$FAO("!6AS !4SL !AS",MODULE,TEXT_LINES,COMMENT)
$ WRITE VERSION VERSION_TEXT
$ SAY VERSION_TEXT
$ SOURCE_LINES = SOURCE_LINES + TEXT_LINES
$ GOTO VERS_FILELOOP
$!
$!   VERSION file is now complete
$!
$VERS_ALLDONE:
$ WRITE VERSION "---------------------------------"
$ WRITE VERSION SOURCE_LINES," lines of source code in ",MODULES," modules"
$ SAY "GAMESS version info is now in SCR:VERSION.LIS"
$ CLOSE VERSION
$ CLOSE DIRLIST
$ DELETE SCR:DIR.LIS;*
$ EXIT

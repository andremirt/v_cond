$! Sep 16, 1990
$!
$ SAY := WRITE SYS$OUTPUT
$!
$ mbldrfile = P1
$ IF mbldrfile.NES."" THEN GOTO GOTNAME
$!
$ SAY " "
$ SAY "Your model builder input is expected to be in file xxx.MB."
$ SAY " "
$ SAY "The file MBLDR.MAN explains what the contents of this input"
$ SAY "ought to be like."
$ SAY " "
$ INQUIRE mbldrfile "Enter xxx (or hit <return> to exit)"
$!
$GOTNAME:
$ mbldrfile = F$EDIT(mbldrfile,"TRIM, UPCASE")
$ IF mbldrfile.EQS."" THEN EXIT
$!
$ ASSIGN/USER_MODE 'mbldrfile'.MB       MBLDR
$ ASSIGN/USER_MODE 'mbldrfile'.LIS      FOR006
$ RUN TOOLS:MBLDR
$!
$ SAY "Your MBLDR output is in file ''mbldrfile'.LIS"
$ EXIT

$! 15 Sep 1990
$!
$ deck = P1
$ IF deck.NES."" THEN GOTO GOTNAME
$!
$ SAY := WRITE SYS$OUTPUT
$ SAY " "
$ SAY "Your CARTIC input should be in file xxx.XYZ."
$ SAY " "
$ SAY "A CARTIC input file contains one card per atom, each card"
$ SAY "with the atom's name, its nuclear charge, and x,y,z coords."
$ SAY " "
$ INQUIRE deck "Enter xxx (or hit <return> to exit)"
$!
$GOTNAME:
$ deck = F$EDIT(deck,"TRIM, UPCASE")
$ IF deck.EQS."" THEN EXIT
$!
$ ASSIGN/USER_MODE 'deck'.XYZ CARTIC
$ ASSIGN/USER_MODE TT:        FOR005
$ ASSIGN/USER_MODE TT:        FOR006
$ RUN TOOLS:CARTIC
$ EXIT

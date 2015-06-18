$! Sep 15, 1990
$!
$ deck = P1
$ IF deck.NES."" THEN GOTO GOTNAME
$!
$ SAY := WRITE SYS$OUTPUT
$ SAY " "
$ SAY "Your old, dirty $VEC group is expected to be in file xxx.VEC."
$ SAY " "
$ SAY "A xxx.VEC file should contain any number of 'comment cards',"
$ SAY "which can be actual GAMESS input cards if you like, followed"
$ SAY "by a ' $VEC' card.  Any preceeding cards will be copied to the"
$ SAY "output file, and will be followed by the clean $VEC group."
$ SAY "Thus your $VEC group should be the LAST THING in xxx.VEC,"
$ SAY "and if it contains fewer MOs than AOs, the input $VEC should"
$ SAY "not contain a $END card.  Be prepared to tell CLENMO how many"
$ SAY "atomic orbitals your molecule is using."
$ SAY " "
$ INQUIRE deck "Enter xxx (or hit <return> to exit)"
$!
$GOTNAME:
$ deck = F$EDIT(deck,"TRIM, UPCASE")
$ IF deck.EQS."" THEN EXIT
$!
$ ASSIGN/USER_MODE 'deck'.VEC        MOIN
$ ASSIGN/USER_MODE SCR:'deck'.MOOUT  MOOUT
$ ASSIGN/USER_MODE TT:               FOR005
$ ASSIGN/USER_MODE TT:               FOR006
$!
$ RUN TOOLS:CLENMO
$ SAY "Your new, cleaned up $VEC is in SCR:''deck'.MOOUT"
$ EXIT

C
C     adapted version for ifort to create actvte.x
C     ifort -o actvte.x actvte.f
C     "*UNX" replaced by 4 spaces
C     andre m 24.1.2012
C
C 25 Jun 04 - MWS - REPLCE taught about lower case substitutions
C 19 MAR 02 - MWS - ADD 4 MORE BLAS ROUTINE NAMES
C 20 DEC 96 - HPP - CHANGES FOR CRAY T3E INTERFACE
C 28 FEB 95 - MWS - FIX COMPLEX PRECISION CHANGES
C 10 JUN 94 - MWS - ADD ONE LEVEL 2 BLAS NAME
C 26 APR 93 - MWS - REMOVE TEST OF MACHINE NAMES
C 20 JUN 91 - MWS - STOP PRINTING THE ACTIVATED LINES
C  7 JUN 89 - MWS - UNIX GETS MACHINE TYPE FROM ENVIRONMENT
C 28 NOV 86 - STE - MAKE BLAS LIST COMPLETE
C 23 OCT 86 - MWS - ACTIVATOR PROGRAM CREATED
C
C     --------------- ACTIVATE GAMESS SOURCE CODE ---------------
C
C     THIS PROGRAM WILL PROCESS GAMESS SOURCE CODE, ACTIVATING
C     THE DESIRED MACHINE DEPENDENT CODE (E.G. THOSE CARDS WITH
C     '*VMS' IN THE FIRST 4 COLUMNS WILL BE CHANGED TO 4 BLANKS.)
C
C     FOR 64 BIT MACHINES, THE REMAINDER OF THE SOURCE LINE
C     IS SCANNED TO CONVERT 'D' EXPONENTS TO 'E' EXPONENTS,
C     TO CHOOSE SINGLE PRECISION FLOATING POINT, AND TO CALL
C     THE SINGLE PRECISION LEVEL 1 BLAS.
C
C     THIS PROGRAM MUST BE --MANUALLY-- ACTIVATED BEFORE COMPILATION.
C     TO RUN ACTVTE ON
C         ALMOST ALL UNIX SYSTEMS,             TURN ON '*UNX' 
C         CRAY T3E,T90,X1 SYSTEM            TURN ON '*UNX' AND '*SNG'
C         NEC-SX                            TURN ON '*UNX' AND '*SNG'
C         VAX OR ALPHA UNDER VMS,              TURN ON '*VMS'
C         IBM MAINFRAME UNDER MVS,             TURN ON '*IBM'
C         IBM MAINFRAME UNDER VM,           DON'T USE THIS PROGRAM.
C     WITHIN ACTVTE ITSELF, BEFORE YOU COMPILE ACTVTE.
C
C     NOTE THAT THE STRING YOU MUST USE TO TURN 'ACTVTE' ON IS
C     NOT NECESSARILY THE SAME STRING WHICH ACTVTE WILL BE 
C     ACTIVATING IN THE GAMESS SOURCE ITSELF.
C
C     FILES USED
C     IN   (INPUT)  - UNACTIVATED GAMESS SOURCE CODE
C     IOUT (OUTPUT) - FORTRAN CODE, READY FOR THE COMPILER
C     IR   (INPUT)  - READ USER SELECTION FOR TARGET MACHINE
C     IW   (OUTPUT) - INFORMATIONAL MESSAGES
C
C     --------------
      PROGRAM ACTVTE
C     --------------
      CHARACTER*4 KEY,BLANK,MACHIN
      CHARACTER*68 TEXT
      CHARACTER*256 FILENM
      LOGICAL FOUND,SINGLE
      DATA BLANK/'    '/
      IN=1
      IOUT=2
      IR=5
      IW=6
C
C     OPEN THE FILES.  IF SYMBOLIC NAMES ARE USED, USE THE NAMES
C     'SRCIN', 'CODEOUT', 'ACTIN', AND 'ACTOUT' FOR ALL MACHINES.
C
      FILENM = ' '
      CALL GETENV('SRCIN',FILENM)
      LENNM=LEN(FILENM)
      IF(LENNM.EQ.0) THEN
         WRITE(IW,*) 'no SRCIN assignment'
         STOP
      END IF
      OPEN(UNIT=IN, FILE=FILENM, STATUS='OLD', ACCESS='SEQUENTIAL',
     *     FORM='FORMATTED')
      FILENM = ' '
      CALL GETENV('CODEOUT',FILENM)
      LENNM=LEN(FILENM)
      IF(LENNM.EQ.0) THEN
         WRITE(IW,*) 'no CODEOUT assignment'
         STOP
      END IF
      OPEN(UNIT=IOUT, FILE=FILENM, STATUS='NEW', ACCESS='SEQUENTIAL',
     *     FORM='FORMATTED')
C
*VMS  OPEN(UNIT=IN,   FILE='SRCIN',   STATUS='OLD', READONLY, SHARED)
*VMS  OPEN(UNIT=IOUT, FILE='CODEOUT', STATUS='NEW')
*VMS  OPEN(UNIT=IR,   FILE='ACTIN',   STATUS='OLD')
*VMS  OPEN(UNIT=IW,   FILE='ACTOUT',  STATUS='NEW')
C
*IBM  OPEN(UNIT=IN,   FILE='SRCIN',   STATUS='OLD')
*IBM  OPEN(UNIT=IOUT, FILE='CODEOUT', STATUS='NEW')
C
C        DETERMINE THE TARGET MACHINE
C
      MACHIN = BLANK
*IBM  MACHIN = '*IBM'
      CALL GETENV('MACHIN',MACHIN)
      IF(MACHIN.EQ.BLANK) READ(IR,900) MACHIN
      WRITE(IW,910) MACHIN
C
      SINGLE=.FALSE.
*SNG  SINGLE=.TRUE.
C
C     LOOP OVER EACH SOURCE CARD
C
      NCARD= 0
      NACT = 0
      NEXP = 0
      NBLAS= 0
      NDBL = 0
  100 CONTINUE
      READ(IN,900,END=800) KEY,TEXT
      NCARD=NCARD+1
C
C     IF KEY IS THE MACHINE REQUESTED, ACTIVATE THIS LINE
C
      IF(KEY.EQ.MACHIN) THEN
C----    WRITE(IW,920) KEY,TEXT
         KEY = BLANK
         NACT=NACT+1
      END IF
      IF(.NOT.SINGLE) GO TO 200
C
C     FOR 64 BIT MACHINES, CONVERT TO SINGLE PRECISION
C
      IF(KEY(1:1).EQ.'C'  .OR.  KEY(1:1).EQ.'*' 
     *                             .AND. KEY(2:4) .NE.'BL3') THEN
      ELSE
         CALL DSCAN(KEY,TEXT,IEXP)
         NEXP = NEXP + IEXP
         CALL TOSNGL(TEXT,NBLAS,NDBL,FOUND)
C----    IF(IEXP.GT.0 .OR. FOUND) WRITE(IW,920) KEY,TEXT
      END IF
C
C     WRITE THE ACTIVATED SOURCE LINE OUT
C
  200 CONTINUE
      WRITE(IOUT,900) KEY,TEXT
      GO TO 100
C
  800 CONTINUE
      WRITE(IW,930) NCARD,NACT
      IF(SINGLE) WRITE(IW,940) NEXP,NBLAS,NDBL
      STOP
C
  900 FORMAT(A4,A68)
  910 FORMAT(1X,'ACTIVATING ',A4,' SOURCE CARDS')
  920 FORMAT(1X,A4,A68)
  930 FORMAT(1X,I6,' CARDS READ,',I4,' CARDS ACTIVATED')
  940 FORMAT(1X,I6,' EXPONENTS WERE CHANGED TO SINGLE PRECISION'/
     *       1X,I6,' BLAS CALLS CHANGED TO SINGLE PRECISION'/
     *       1X,I6,' DOUBLE PRECISIONS CHANGED TO TYPE REAL')
      END
C     -------------------------------
      SUBROUTINE DSCAN(KEY,TEXT,IEXP)
C     -------------------------------
      CHARACTER*(*) KEY,TEXT
      CHARACTER*1   PLUS,MINUS,DOT,STAR,LETC,KHAR
      CHARACTER*10  NUMS
      DATA PLUS/'+'/, MINUS/'-'/, DOT/'.'/, STAR/'*'/, LETC/'C'/
      DATA NUMS/'0123456789'/
C
C     THIS ROUTINE EXAMINES ALL OCCURENCES OF THE LETTER 'D'
C     IN THE SOURCE LINE, TO SEE IF THE 'D' IS THE EXPONENT
C     PORTION OF A FLOATING POINT CONSTANT.  IF THE 'D' IS
C     FOLLOWED BY A SIGN, AND PRECEEDED BY NUMBERS AND A
C     DECIMAL POINT, THE 'D' IS CONVERTED TO A 'E'.
C
      IEXP=0
C
C     PASS COMMENT LINES STRAIGHT THROUGH
C
      KHAR = KEY(1:1)
      IF(KHAR.EQ.LETC) RETURN
      IF(KHAR.EQ.STAR) RETURN
C
C     LOOP OVER ALL 'D' CHARACTERS IN THE LINE
C
      KOL=0
  100 CONTINUE
      KOL=KOL+1
      CALL NEXTD(TEXT,KOL)
      IF(KOL.EQ.0) RETURN
C
C     CHECK TO THE RIGHT, 'D' MUST BE FOLLOWED BY A SIGNED
C     EXPONENT FOR THIS TO ACTUALLY BE AN EXPONENT.
C
      KK = KOL+1
      KHAR = TEXT(KK:KK)
      IF(KHAR.EQ. PLUS) GO TO 200
      IF(KHAR.EQ.MINUS) GO TO 200
      GO TO 100
C
  200 CONTINUE
      KK = KK+1
      KHAR = TEXT(KK:KK)
      IF(INDEX(NUMS,KHAR).EQ.0) GO TO 100
C
C     CHECK TO THE LEFT, DECIMAL POINT IS ACCEPTABLE,
C     NUMBER FOUND CONTINUES LEFTWARD CHECK, ANYTHING
C     ELSE MEANS THIS 'D' IS NOT AN EXPONENT.
C
      KK = KOL
  300 CONTINUE
         KK = KK-1
         IF(KK.LT.1) GO TO 100
         KHAR = TEXT(KK:KK)
         IF(KHAR.EQ.DOT) GO TO 400
         IF(INDEX(NUMS,KHAR).EQ.0) GO TO 100
         GO TO 300
C
C     VALID EXPONENT FOUND, CHANGE IT TO 'E'
C
  400 CONTINUE
      TEXT(KOL:KOL)='E'
      IEXP = IEXP+1
      GO TO 100
      END
C     ---------------------------
      SUBROUTINE NEXTD(TEXT,ILOC)
C     ---------------------------
      CHARACTER TEXT*(*)
C
C     LOOK FOR NEXT LETTER D, STARTING FROM COLUMN ILOC
C
      L1=LEN(TEXT)
      JLOC=INDEX(TEXT(ILOC:L1),'D')
      IF(JLOC.GT.0) THEN
            ILOC=ILOC+JLOC-1
         ELSE
            ILOC=0
         ENDIF
      RETURN
      END
C     ----------------------------------------
      SUBROUTINE TOSNGL(TEXT,NBLAS,NDBL,FOUND)
C     ----------------------------------------
      CHARACTER*65 TEXT
      LOGICAL FOUND,FNDBLA,FNDDBL
C
C         CHANGE DOUBLE PRECISION NAMES FOR BLAS ROUTINES
C         TO THEIR SINGLE PRECISION COUNTERPARTS.
C
      FNDBLA=.FALSE.
C         COMPLETE LIST OF LEVEL 1 BLAS
      CALL REPLCE(TEXT,'DASUM' ,'SASUM' ,FNDBLA)
      CALL REPLCE(TEXT,'DAXPY' ,'SAXPY' ,FNDBLA)
      CALL REPLCE(TEXT,'DCOPY' ,'SCOPY' ,FNDBLA)
      CALL REPLCE(TEXT,'DDOT'  ,'SDOT'  ,FNDBLA)
      CALL REPLCE(TEXT,'DNRM2' ,'SNRM2' ,FNDBLA)
      CALL REPLCE(TEXT,'DROT'  ,'SROT'  ,FNDBLA)
      CALL REPLCE(TEXT,'DROTG' ,'SROTG' ,FNDBLA)
      CALL REPLCE(TEXT,'DSCAL' ,'SSCAL' ,FNDBLA)
      CALL REPLCE(TEXT,'DSWAP' ,'SSWAP' ,FNDBLA)
      CALL REPLCE(TEXT,'IDAMAX','ISAMAX',FNDBLA)
      CALL REPLCE(TEXT,'ZDOTC' ,'CDOTC' ,FNDBLA)
      CALL REPLCE(TEXT,'ZSWAP' ,'CSWAP' ,FNDBLA)
      CALL REPLCE(TEXT,'DZNRM2','SCNRM2',FNDBLA)
C         THREE LEVEL 2 ROUTINES ARE USED
      CALL REPLCE(TEXT,'DGER'  ,'SGER'  ,FNDBLA)
      CALL REPLCE(TEXT,'DTRMV' ,'STRMV' ,FNDBLA)
      CALL REPLCE(TEXT,'DGEMV' ,'SGEMV' ,FNDBLA)
C         FOUR LEVEL 3 ROUTINES ARE USED
      CALL REPLCE(TEXT,'DTRMM' ,'STRMM' ,FNDBLA)
      CALL REPLCE(TEXT,'DTRSM' ,'STRSM' ,FNDBLA)
      CALL REPLCE(TEXT,'DGEMM' ,'SGEMM' ,FNDBLA)
      CALL REPLCE(TEXT,'ZGEMM' ,'CGEMM' ,FNDBLA)
C         LAPACK SUBSTITUTION
      CALL REPLCE(TEXT,'ZHEEV' ,'CHEEV' ,FNDBLA)
      IF(FNDBLA) NBLAS=NBLAS+1
C
      FNDDBL=.FALSE.
      CALL REPLCE(TEXT,'DOUBLE PRECISION',
     *                 'REAL            ',FNDDBL)
C         THERE ARE SOME COMPLEX NUMBERS IN SPIN-ORBIT
      CALL REPLCE(TEXT,'COMPLEX*16','COMPLEX   ',FNDDBL)
      CALL REPLCE(TEXT,'DCMPLX'    ,' CMPLX'    ,FNDDBL)
      CALL REPLCE(TEXT,'DREAL'     ,' REAL'     ,FNDDBL)
      CALL REPLCE(TEXT,'DIMAG'     ,'AIMAG'     ,FNDDBL)
      CALL REPLCE(TEXT,'DCONJG'    ,' CONJG'    ,FNDDBL)
      IF(FNDDBL) NDBL=NDBL+1
C
      FOUND = FNDBLA.OR.FNDDBL
      RETURN
      END
C     -------------------------------------------
      SUBROUTINE REPLCE(TEXT,STRNG1,STRNG2,FOUND)
C     -------------------------------------------
      CHARACTER STRNG1*(*),STRNG2*(*)
      CHARACTER*65 TEXT,strng
      LOGICAL FOUND
C
C           LOCATE ALL OCCURENCES OF 'STRNG1' IN THE
C           STRING 'TEXT', CHANGE THEM TO 'STRNG2'.
C       'STRNG1' AND 'STRNG2' MUST HAVE THE SAME LENGTH
C
      LEN1=LEN(STRNG1)
      LEN2=LEN(STRNG2)
      IF(LEN2.NE.LEN1) THEN
         WRITE(6,*) 'ILLEGAL STRING CHANGE REQUESTED'
         STOP
      END IF
c
c        ipass=1 looks for caps, which is how STRNG1 should be given.
c        ipass=2 looks for lower case.
c
      do 200 ipass=1,2
         if(ipass.eq.1) then
            strng(1:len1) = strng1(1:len1)
         else
            do k=1,len1
               num = ichar(strng1(k:k))
               if(num.ge.65 .and. num.le.90) num=num+32
               strng(k:k)=char(num)
            enddo
         end if
c
c    the string might occur more than once, so loop until not found
c
  100    CONTINUE
            LOC = INDEX(TEXT,STRNG(1:len1))
            IF(LOC.EQ.0) go to 200
            TEXT(LOC:LOC+LEN1-1) = STRNG2
            FOUND = .TRUE.
         GO TO 100
  200 continue
      return
C
C          some but not all Cray models like the Posix interface
C
*PSX  END
*PSX  SUBROUTINE GETENV(NAME,VALUE)
*PSX  CHARACTER*(*) NAME,VALUE
*PSX  LENNAM = LEN(NAME)
*PSX  CALL PXFGETENV(NAME,LENNAM,VALUE,LENVAL,IERROR)
*PSX  RETURN
      END

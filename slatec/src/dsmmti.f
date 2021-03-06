*DECK DSMMTI
      SUBROUTINE DSMMTI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***BEGIN PROLOGUE  DSMMTI
C***PURPOSE  SLAP MSOLVE for LDU Factorization of Normal Equations.
C            This routine acts as an interface between the SLAP generic
C            MMTSLV calling convention and the routine that actually
C                                    -1
C            computes  [(LDU)*(LDU)']  B = X.
C***LIBRARY   SLATEC (SLAP)
C***CATEGORY  D2E
C***TYPE      DOUBLE PRECISION (SSMMTI-S, DSMMTI-D)
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
C***AUTHOR  Greenbaum, Anne, (Courant Institute)
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-60
C             Livermore, CA 94550 (510) 423-3141
C             seager@llnl.gov
C***DESCRIPTION
C       It is assumed that RWORK and IWORK have initialized with
C       the information required for DSMMI2:
C          IWORK(1) = Starting location of IL in IWORK.
C          IWORK(2) = Starting location of JL in IWORK.
C          IWORK(3) = Starting location of IU in IWORK.
C          IWORK(4) = Starting location of JU in IWORK.
C          IWORK(5) = Starting location of L in RWORK.
C          IWORK(6) = Starting location of DINV in RWORK.
C          IWORK(7) = Starting location of U in RWORK.
C       See the DESCRIPTION of DSMMI2 for details.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DSMMI2
C***REVISION HISTORY  (YYMMDD)
C   871119  DATE WRITTEN
C   881213  Previous REVISION DATE
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
C   890922  Numerous changes to prologue to make closer to SLATEC
C           standard.  (FNF)
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)
C   910411  Prologue converted to Version 4.0 format.  (BAB)
C   920511  Added complete declaration section.  (WRB)
C   921113  Corrected C***CATEGORY line.  (FNF)
C   930701  Updated CATEGORY section.  (FNF, WRB)
C***END PROLOGUE  DSMMTI
C     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
C     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(10), JA(NELT)
C     .. Local Scalars ..
      INTEGER LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU
C     .. External Subroutines ..
      EXTERNAL DSMMI2
C***FIRST EXECUTABLE STATEMENT  DSMMTI
C
C         Pull out the locations of the arrays holding the ILU
C         factorization.
C
      LOCIL = IWORK(1)
      LOCJL = IWORK(2)
      LOCIU = IWORK(3)
      LOCJU = IWORK(4)
      LOCL = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU = IWORK(7)
C
      CALL DSMMI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL),
     $     RWORK(LOCL), RWORK(LOCDIN), IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU))
C
      RETURN
C------------- LAST LINE OF DSMMTI FOLLOWS ----------------------------
      END

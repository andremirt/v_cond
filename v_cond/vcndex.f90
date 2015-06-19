program vcndex
  !***********************************************************************
  !
  ! Program to calculate values of conditional potential in a number
  !  of points. The points usually represent a numerical integration
  !  grid and are read from file "points". The calculation of vcond in all
  !  points (x,y,z) involves integrals of the type :  phi(i)*phi(j)/r12,
  !  integrate over 2. These integrals are of the same type as the
  !  electron-nuclear attraction integrals and the phi's should be 
  !  cartesian gaussians.
  !
  ! the program can handle Gamess-US densities and basis sets
  ! it has to be linked with the SLATEC program library, see below
  !
  !***********************************************************************
  use integral
  use readgauss
  implicit none

  integer           :: norb, idmp, ip2un, allst
  double precision  :: rnel
  real              :: starttime, endtime
  character(8)      :: date
  character(10)     :: time

  open( 6, file='vcond.output')
  write(6,'(''***********************************************************************************'')')
  write(6,'(''***                                                                             ***'')')
  write(6,'(''***  This is program vcndex to calculate the conditional and Hartree potential  ***'')')
  write(6,'(''***                                                                             ***'')')
  write(6,'(''***       date of program build   :   unknowndate             ***'')')
  write(6,'(''***                                                                             ***'')')
  write(6,'(''***********************************************************************************'')')

  call date_and_time( date, time)
  write(6,'(/''       begining with program execution on    '',a2,1x,a2,1x,a4,3x,a2,'':'',a2,'':'',a2/)') &
       & date(5:6), date(7:8), date(1:4), time(1:2), time(3:4), time(5:6)
  write(6,'(''***********************************************************************************'')')
  call cpu_time( starttime)

  !unit numbers of files
  idmp = 8
  ip2un = 7

  norb = 0
  rnel = 0.d0

  open( idmp, file='vcond.basis', status='old', err=223)
  rewind(idmp)

  open( ip2un, file='vcond.input', status='old', err=223)
  rewind( ip2un)
  read(ip2un,*,err=224) rnel, norb

  if (norb.le.0) then
     write(6,'(/'' ERROR; norb not properly specified'')')
     stop
  endif
  write(6,'(/''  norb  = '',i4)') norb

  if (rnel.le.0.d0) then
     write(6,'(/'' ERROR; rnel not properly specified'')')
     stop
  endif
  write(6,'(''  # of electrons : '',f8.2)') rnel

  call getgss( idmp)
  call cOlAO
  
  call vbrain( norb, ip2un, idmp, rnel)

  call date_and_time( date, time)
  write(6,'(/''***********************************************************************************'')')
  write(6,'(/''       finishing with program execution on    '',a2,1x,a2,1x,a4,3x,a2,'':'',a2,'':'',a2)') &
       & date(5:6), date(7:8), date(1:4), time(1:2), time(3:4), time(5:6)
  call cpu_time( endtime)
  write(6,'(''       elapsed time in sec : '',f12.2/)') endtime - starttime
  write(6,'(''***********************************************************************************'')')
  write(6,'(''***********************************************************************************'')')

  close(6)
  close(ip2un)
  close(idmp)

  stop

223 write(6,'(/'' ERROR; file "input" does not exist'')')
  close(6)
  close(idmp)
  stop

224 write(6,'(/'' ERROR; number of electrons or orbitals not properly given'')')
  close(6)
  stop

end program vcndex
!!$**********************************************************************************
!!$
!!$ this is the function d1mach for the slatec library writen in Fortran90 syntax
!!$ you might need it in case of trouble with compiling slatec
!!$
!!$**********************************************************************************
!!$REAL(KIND=8) FUNCTION D1MACH(I)
!!$
!!$  INTEGER,INTENT(IN):: I
!!$  !
!!$  !  DOUBLE-PRECISION MACHINE CONSTANTS
!!$  !  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!!$  !  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!!$  !  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!!$  !  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!!$  !  D1MACH( 5) = LOG10(B)
!!$  !                    
!!$  REAL(KIND=8),DIMENSION(5),SAVE::X
!!$  LOGICAL,SAVE:: ran = .false.
!!$
!!$  IF(.NOT. ran)THEN
!!$
!!$     X(1) = TINY(X(1))
!!$     X(2) = HUGE(X(2))
!!$     X(3) = EPSILON(X(3))/RADIX(X(3))
!!$     X(4) = EPSILON(X(4))
!!$     X(5) = LOG10(REAL(RADIX(X(5)),KIND=8))
!!$     ran = .true.
!!$
!!$  END IF
!!$
!!$  IF (I < 1 .OR. I > 5) THEN
!!$
!!$     WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
!!$     STOP
!!$
!!$  END IF
!!$
!!$  D1MACH = X(I)
!!$
!!$  RETURN
!!$
!!$END FUNCTION D1MACH
!!$
!!$*****************************************************************************
!!$
!!$ this is the botom of the makefile for the static libraries
!!$
!!$*****************************************************************************
!!$ #
!!$ #  Static library target
!!$ #
!!$
!!$ libslatec.a:  $(OBJ) d1mach
!!$ 	ar r $@ $(OBJ) d1mach.o; ranlib $@
!!$
!!$ ${OBJ}: %.o:
!!$ 	$(FC) $(FFLAGS) -c ../src/$(@:.o=.f) -o $@
!!$
!!$ d1mach: d1mach.o
!!$ 	${FC} ${FFLAGS} -c d1mach.f90 -o d1mach.o
!!$
!!$ clean:
!!$ 	rm -f *.o libslatec.a


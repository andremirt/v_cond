module readgauss

!  use integral
  implicit none
  private
  public bset, geom, getgss, nat, natt, AngMTab, Pi, nbsf, cgauss, tgeom, tAngMTab, pgauss
    
  type pgauss
     integer                       :: nprims, l
     double precision, allocatable :: coeff(:), coeff2(:), exp(:), NormCoeff(:,:)
  end type pgauss
  
  type cgauss
     character                 :: name*10
     integer                   :: nshells
     type(pgauss), allocatable :: primtvs(:)
  end type cgauss
  
  type(cgauss), allocatable, save :: bset(:)
  
  type tAngMTab
     integer                       :: ncomb
     integer, allocatable          :: lmn(:,:)
  end type tAngMTab

  type(tAngMTab), save :: AngMTab(6)

  type tgeom
     character        :: name*10
     double precision :: coord(3)
  end type tgeom
  
  type(tgeom), allocatable, save :: geom(:)

  integer,save                :: nat, natt, nbsf
  integer                     :: lmax, allst
  double precision, external  :: dfac
  double precision, parameter :: Pi = 3.14159265358979d0

contains
  
  subroutine getgss(idmp)
    !*************************************************************
    !
    ! here we read the information to the contracted gaussians
    ! gaussians should be cartesian and normalized
    !
    !*************************************************************
    integer, intent(in) :: idmp

    integer :: i, j, k
    
    write(6,'(/'' reading geometry and basis set input from file "basis"'')')    

    AngMTab(1)%ncomb = 1
    AngMTab(2)%ncomb = 3
    AngMTab(3)%ncomb = 6
    AngMTab(4)%ncomb = 10
    AngMtab(5)%ncomb = 15
    AngMtab(6)%ncomb = 1

    lmax = 0
    nbsf = 0

    !read geometry from input in cartesian form 
    read(idmp,*) nat
    write(6,'(''  atoms found in basis file : '',i4)') nat

    allocate( geom(nat), stat = allst)
    if( allst .ne. 0) then
       write(6,'(/''ERROR in allocation: geom!!!'')')
       stop
    end if

    write(6,'(''  the geometry is :''/)')

    do i =  1, nat
       read(idmp,*) geom(i)%name, geom(i)%coord
       write(6,'(4x,a10,3f12.4)') geom(i)%name, geom(i)%coord
    end do

    !read basisset information from input
    read(idmp,*) natt

    if ( nat .ne. natt) then
       write(6,'(/''ERROR: number of atoms does not match to number of basissets given!!!'')')
       stop
    end if
    
    allocate( bset(natt), stat = allst)
    if( allst .ne. 0) then
       write(6,'(/''ERROR in allocation: bset!!!'')')
       stop
    end if

    write(6,'(/''  the basis set is :''/)')

    do i = 1, natt

       read(idmp,*) bset(i)%name
       write(6,'(4x,a10)') bset(i)%name

       if ( bset(i)%name .ne. geom(i)%name) then
          write(6,'(/''ERROR: basisset for atom'',i5,'' not properly given!!!'')') i
          write(6,'(/''expected set :'',a12,'', given set :'',a12)') geom(i)%name, bset(i)%name
          stop
       end if

       read(idmp,*) bset(i)%nshells
       allocate( bset(i)%primtvs(bset(i)%nshells))

       do j = 1, bset(i)%nshells

          read(idmp,*) bset(i)%primtvs(j)%l, bset(i)%primtvs(j)%nprims
          nbsf = nbsf + AngMTab(bset(i)%primtvs(j)%l)%ncomb
          write(6,'(6x,''shell: '',i4,'' angular momentum: '',i4)') j, bset(i)%primtvs(j)%l - 1

!!$          if ( bset(i)%primtvs(j)%l .eq. 8) then
!!$             !this is the case for sp-shells
!!$
!!$             if ( lmax .lt. 2) lmax = 2
!!$
!!$             ntshells = ntshells + 3
!!$             allocate( bset(i)%primtvs(j)%coeff(bset(i)%primtvs(j)%nprims),&
!!$                  & bset(i)%primtvs(j)%exp(bset(i)%primtvs(j)%nprims),&
!!$                  & bset(i)%primtvs(j)%coeff2(bset(i)%primtvs(j)%nprims))
!!$             
!!$             do k = 1, bset(i)%primtvs(j)%nprims
!!$
!!$                read(idmp,*) bset(i)%primtvs(j)%exp(k), &
!!$                     & bset(i)%primtvs(j)%coeff(k), bset(i)%primtvs(j)%coeff2(k)
!!$                
!!$             enddo
!!$             
!!$             cycle
!!$
!!$          endif

          if ( bset(i)%primtvs(j)%l .gt. lmax) lmax = bset(i)%primtvs(j)%l

          allocate( bset(i)%primtvs(j)%coeff(bset(i)%primtvs(j)%nprims),&
               & bset(i)%primtvs(j)%exp(bset(i)%primtvs(j)%nprims), &
               & bset(i)%primtvs(j)% &
               & NormCoeff(bset(i)%primtvs(j)%nprims,AngMTab(bset(i)%primtvs(j)%l)%ncomb))

          do k = 1, bset(i)%primtvs(j)%nprims

             read(idmp,*) bset(i)%primtvs(j)%exp(k), &
                  & bset(i)%primtvs(j)%coeff(k)
             write(6,'(6x,f12.6,f14.8)') bset(i)%primtvs(j)%exp(k), &
                  & bset(i)%primtvs(j)%coeff(k)

          enddo

       enddo

    enddo

    write(6,*)
    write(6,'(''  total number of AOs : '',i5)') nbsf
    write(6,*)

    call cAngMTab( lmax)
    call normzAO

    return
    
  end subroutine getgss

  subroutine cAngMTab( lmax)
    !**********************************************************************
    !
    ! here for each angular momentum L the powers lmn of the cartesian 
    !  gausians are defined
    !
    !
    !      (x-Ax)**lmn(L,1) * (y-Ay)**lmn(L,2) * (z-Az)**lmn(L,3)
    !
    !
    ! the ordering of the powers corresponds to the ordering 
    !  of Gamess-US
    !
    !**********************************************************************
    integer, intent(in) :: lmax

    integer :: i

    if ( lmax .gt. 4) then
       write(6,'(''ERROR: support only up to f-type basis functions'')')
       write(6,'(''       see subroutine cAngMTab'')')
       stop
    end if

    allocate( AngMTab(1)%lmn(3,1), AngMTab(6)%lmn(3,4), stat = allst)
    if( allst .ne. 0) then
       write(6,'(''ERROR in allocation: AngMTab%lmn'')')
       stop
    end if

    !s-shell
    AngMTab(1)%lmn(:,1) = 0
    
    if ( lmax .eq. 1) return
    
    allocate( AngMTab(2)%lmn(3,3), stat = allst)
    if( allst .ne. 0) then
       write(6,'(''ERROR in allocation: AngMTab%lmn'')')
       stop
    end if

    !p-shell
    AngMTab(2)%lmn(:,:) = 0
    AngMTab(2)%lmn(1,1) = 1
    AngMTab(2)%lmn(2,2) = 1
    AngMTab(2)%lmn(3,3) = 1

    !sp-shell
    AngMTab(6)%lmn(:,:) = 0
    AngMTab(6)%lmn(1,2) = 1
    AngMTab(6)%lmn(2,3) = 1
    AngMTab(6)%lmn(3,4) = 1
    
    if ( lmax .eq. 2) return
    
    allocate( AngMTab(3)%lmn(3,6), stat = allst)
    if( allst .ne. 0) then
       write(6,'(''ERROR in allocation: AngMTab%lmn'')')
       stop
    end if

    !d-shell
    AngMTab(3)%lmn(:,:) = 0
    AngMTab(3)%lmn(1,1) = 2
    AngMTab(3)%lmn(2,2) = 2
    AngMTab(3)%lmn(3,3) = 2
    AngMTab(3)%lmn(1,4) = 1
    AngMTab(3)%lmn(2,4) = 1
    AngMTab(3)%lmn(1,5) = 1
    AngMTab(3)%lmn(3,5) = 1
    AngMTab(3)%lmn(2,6) = 1
    AngMTab(3)%lmn(3,6) = 1

    if ( lmax .eq. 3) return
    
    allocate( AngMTab(4)%lmn(3,10), stat = allst)
    if( allst .ne. 0) then
       write(6,'(''ERROR in allocation: AngMTab%lmn'')')
       stop
    end if

    !f-shell
    AngMTab(4)%lmn(:,:) = 0
    AngMTab(4)%lmn(1,1) = 3
    AngMTab(4)%lmn(2,2) = 3
    AngMTab(4)%lmn(3,3) = 3
    AngMTab(4)%lmn(1,4) = 2
    AngMTab(4)%lmn(2,4) = 1
    AngMTab(4)%lmn(1,5) = 2
    AngMTab(4)%lmn(3,5) = 1
    AngMTab(4)%lmn(2,6) = 2
    AngMTab(4)%lmn(1,6) = 1
    AngMTab(4)%lmn(2,7) = 2
    AngMTab(4)%lmn(3,7) = 1
    AngMTab(4)%lmn(3,8) = 2
    AngMTab(4)%lmn(1,8) = 1
    AngMTab(4)%lmn(3,9) = 2
    AngMTab(4)%lmn(2,9) = 1
    AngMTab(4)%lmn(1,10) = 1
    AngMTab(4)%lmn(2,10) = 1
    AngMTab(4)%lmn(3,10) = 1

    if ( lmax .eq. 4) return
    
    allocate( AngMTab(5)%lmn(3,15), stat = allst)
    if( allst .ne. 0) then
       write(6,'(''ERROR in allocation: AngMTab%lmn'')')
       stop
    end if

    !g-shell
    AngMTab(5)%lmn(:,:) = 0
    AngMTab(5)%lmn(1,1) = 4
    AngMTab(5)%lmn(2,2) = 4
    AngMTab(5)%lmn(3,3) = 4
    AngMTab(5)%lmn(1,4) = 3
    AngMTab(5)%lmn(2,4) = 1
    AngMTab(5)%lmn(1,5) = 3
    AngMTab(5)%lmn(3,5) = 1
    AngMTab(5)%lmn(2,6) = 3
    AngMTab(5)%lmn(1,6) = 1
    AngMTab(5)%lmn(2,7) = 3
    AngMTab(5)%lmn(3,7) = 1
    AngMTab(5)%lmn(3,8) = 3
    AngMTab(5)%lmn(1,8) = 1
    AngMTab(5)%lmn(3,9) = 3
    AngMTab(5)%lmn(2,9) = 1
    AngMTab(5)%lmn(1,10) = 2
    AngMTab(5)%lmn(2,10) = 2
    AngMTab(5)%lmn(1,11) = 2
    AngMTab(5)%lmn(3,11) = 2
    AngMTab(5)%lmn(2,12) = 2
    AngMTab(5)%lmn(3,12) = 2
    AngMTab(5)%lmn(1,13) = 2
    AngMTab(5)%lmn(2,13) = 1
    AngMTab(5)%lmn(3,13) = 1
    AngMTab(5)%lmn(1,14) = 1
    AngMTab(5)%lmn(2,14) = 2
    AngMTab(5)%lmn(3,14) = 1
    AngMTab(5)%lmn(1,15) = 1
    AngMTab(5)%lmn(2,15) = 1
    AngMTab(5)%lmn(3,15) = 2

    return

  end subroutine cAngMTab

  subroutine normzAO

    integer          :: iat, ish, ipr, iam, lmn(3)
    double precision :: alpha

    write(6,'('' normalizing cartesian gaussians'')')

    do iat = 1, nat

       do ish = 1, bset(iat)%nshells

          do iam = 1, AngMTab(bset(iat)%primtvs(ish)%l)%ncomb

             lmn = AngMTab(bset(iat)%primtvs(ish)%l)%lmn(:,iam)

             do ipr = 1, bset(iat)%primtvs(ish)%nprims

                alpha = bset(iat)%primtvs(ish)%exp(ipr)
                bset(iat)%primtvs(ish)%NormCoeff(ipr,iam) =  ( 2 * alpha / Pi)**( 3.d0 / 4) * &
                     & dsqrt( ( 8 * alpha)**( sum( lmn)) * dfac( lmn(1)) * dfac( lmn(2)) * &
                     & dfac( lmn(3)) / ( dfac( 2 * lmn(1)) * dfac( 2 * lmn(2)) &
                     & * dfac( 2 * lmn(3)))) * bset(iat)%primtvs(ish)%coeff(ipr)

             end do
          end do
       end do
    end do
    
    return

  end subroutine normzAO
  
end module readgauss

module density

  implicit none
  public ijkl, pdfmo, p1dens, pmo

  integer, allocatable, save          :: ijkl(:,:)
  double precision, allocatable, save :: pdfmo(:), pmo(:)

contains

  subroutine p1dens( iunit, norb, rnel, nbsf)
    !*******************************************************************
    !
    ! reads pair-distribution function in MO basis and
    !  calculates density in MO basis
    !
    ! density stored in array pmo in triangular form (that is: the sum of
    !  of the off-diagonal elements (i,j) and (j,i) is stored in position
    !  pmo(i*(i-1)/2+j),  assuming i>j)
    !
    ! arguments:
    !  iunit : unit-number of two-matrix file
    !  norb  : number of orbitals on two-matrix file 
    !  rnel  : number of electrons
    !
    !*******************************************************************
    integer, intent(in)           :: iunit, norb, nbsf
    double precision, intent(in)  :: rnel

    integer          :: imb, nmb, ntot, ntot2, ialloc, m, i, j, k, l, &
         & ij, kl, imo, iimo
    double precision :: normpdf, normp, fctr

    write(6,'(/'' reading pdf in MO basis from file "input"''/'' and calculating density in MO basis'')')

    imb = 0
    nmb = 0
    normp = 0.d0
    normpdf = 0.d0
    ntot = norb * ( norb + 1) / 2
    ntot2 = ntot * ( ntot + 1) / 2

    allocate( pmo(nbsf*(nbsf+1)/2), pdfmo(ntot2), ijkl(4,ntot2), stat = ialloc)
    if ( ialloc .ne. 0) then 
       write(6,'(/'' No memory allocation pmo, pdfmo'')')
       stop
    endif

    pmo = 0.d0
    pdfmo = 0.d0
    ijkl = 0

    !read dm2 in vector p2mo 
    read(iunit,*)

    do i = 1, ntot2
 
       imb = imb + 1
       read(iunit,'(4i5,es24.16)',err=10) ijkl( 1:4, i), pdfmo(i)

    enddo
    !
    !**   normalize dm2 to N(N-1)/2 instead of N(N-1)
    !
    !      p2mo = p2mo / 2
    !
    !*** loop over p2mo elements in this block
    !
    do m = 1, ntot2

       nmb = nmb + 1
       i = ijkl(1,m)
       j = ijkl(2,m)
       k = ijkl(3,m)
       l = ijkl(4,m)

       if ( (i.lt.1) .or. (j.lt.1) .or. (k.lt.1) .or. (l.lt.1) .or. (i.gt.norb) &
            & .or. (j.gt.norb) .or. (k.gt.norb) .or. (l.gt.norb)) then

          write(6,'(/,'' ERROR; nonexistent determinant '',i3,1x,i3,1x,i3,1x,i3, &
               & '' encountered'')') i, j, k, l
          stop

       else

          ij = i * ( i - 1) / 2 + j
          kl = k * ( k - 1) / 2 + l
!          if ( ij .eq. kl) write(*,*) 'ij eq kl', ijkl(:,m)

          !calculates density from pdf
          !(transformation from pdf matrix of Gamess routine MATRD2 format to unfolded format,
          !but it's not done here since internal format is different)
          if ( i .eq. j .and. k .eq. l) normpdf = normpdf + pdfmo(m)

          if ( i .eq. j .and. i .eq. k .and. i .eq. l) then
             pmo(ij) = pmo(ij) + pdfmo(m)
!             pdfmo(m) = pdfmo(m)
          elseif ( j .eq. k .and. j .eq. l) then
!             write(*,*) 'ijjj', ijkl(:,m)
             pmo(ij) = pmo(ij) + pdfmo(m) * .5d0
!             pdfmo(m) = pdfmo(m) / 4
          elseif ( i .eq. k .and. j .eq. l) then
!             pdfmo(m) = pdfmo(m) / 4
          elseif ( i .eq. j .and. k .eq. l) then
             pmo(ij) = pmo(ij) + pdfmo(m) * .5d0
             pmo(kl) = pmo(kl) + pdfmo(m) * .5d0
!             pdfmo(m) = pdfmo(m) / 2
          elseif ( i .eq. j .and. i .eq. k) then
             pmo(kl) = pmo(kl) + pdfmo(m) * .5d0
!             pdfmo(m) = pdfmo(m) / 4
          elseif ( j .eq. l) then
!             pdfmo(m) = pdfmo(m) / 8
          elseif ( k .eq. l) then
             pmo(ij) = pmo(ij) + pdfmo(m) * .5d0
!             pdfmo(m) = pdfmo(m) / 4
          elseif ( j .eq. k) then
!             pdfmo(m) = pdfmo(m) / 8
          elseif ( i .eq. k) then
!             pdfmo(m) = pdfmo(m) / 8
          elseif ( i .eq. j) then
             pmo(kl) = pmo(kl) + pdfmo(m) * .5d0
!             pdfmo(m) = pdfmo(m) / 4
          else
!             pdfmo(m) = pdfmo(m) / 8
          end if

       endif

    enddo
    
    !transformation to internal format as used in subroutine cond
    pdfmo = pdfmo * .5d0
    
!!$    do m = 1, ntot2
!!$
!!$       write(*,*) 'pdfmo', ijkl(:,m), pdfmo(m)
!!$
!!$    enddo
!!$
!!$    do m = 1, ntot
!!$
!!$       write(*,*) 'pmo', pmo(m)
!!$
!!$    enddo
    !
    !***  normalize density to N
    !     assumption: pdf normalized to N(N-1)
    !
    fctr = 1.0d0 / ( rnel - 1.d0)
    pmo = pmo * fctr

    do imo = 1, norb

       iimo = imo * ( imo + 1) / 2
       normp = normp + pmo(iimo)

    enddo

    write(6,'(''    ====> complete'')')
    write(6,'(''  number of p2(ijkl)-elements restored from file : '', &
         & i12,2x,i12)') nmb, imb

    if ( abs( normpdf - rnel*(rnel-1)) .gt. 1.d-6) then
       write(6,'(''  WARNING!!! pdf not proper normalized to N(N-1)'')')
       write(6,'(''  N(N- 1) nominal: '',f16.10,'',  actual : '',g17.10)') &
            & rnel * ( rnel - 1), normpdf
       write(6,'(''  diff : '',g17.6)') abs( normpdf - rnel * ( rnel - 1))
    endif

    write(6,'(''  pdf is normalized to N(N-1) : '',f14.8)') normpdf

    if ( abs( normp - rnel) .gt. 1.d-6) then
       write(6,'(''  WARNING!!! density not proper normalized to N'')')
       write(6,'(''  N nominal: '',f16.10,'',  actual: '',g17.10)') rnel, normp
       write(6,'(''  diff : '',g17.6)') abs( rnel - normp)
    endif

    write(6,'(''  density is normalized to  N : '',f14.8)') normp

    return

10  continue
    write(6,'(/,'' ERROR::: something went wrong reading pdf'')')
    write(6,'(''pdf elements expected : '',i8,/,''pdf element attempted to read : '',i8)') ntot2, imb
    stop

  end subroutine p1dens

end module density

subroutine vbrain( norb, ip2un, idmp, rnel)

  use integral
  use readgauss
  use density
  use mvalAO

  implicit none

  integer, intent(in)          :: norb, ip2un, idmp
  double precision, intent(in) :: rnel

  integer          :: i, ii, j, ij, npnt, ipnt
  double precision :: pao(nbsf*(nbsf+1)/2), valAO(nbsf), valMO(nbsf), &
       & vmoao(nbsf,nbsf), Tvmoao(nbsf,nbsf), Pcoord(3), vhartr, vcond, rhovcond, &
       & rpcond(nbsf*(nbsf+1)/2), rpcndao(nbsf*(nbsf+1)/2), densty, &
       & NATAOInt(nbsf*(nbsf+1)/2)
!       & NATAOInt(nbsf,nbsf)
!       & XAOInt(nbsf,nbsf), rhovX, vX

  !initialise files, and read number of points
  open( 11, file = 'vhartr.dat')
  open( 12, file = 'vcond.dat')
  open( 13, file = 'density.dat')
  open( 14, file = 'rhovcond.dat')
  open( 15, file = 'wXC.dat')
!  open( 16, file = 'vX.dat')
  rewind(11)
  rewind(12)
  rewind(13)
  rewind(14)
  rewind(15)
!  rewind(16)

  !read pair-distribution function in MO basis and 
  !calculate density in MO basis from pdf
  call p1dens( ip2un, norb, rnel, nbsf)
  !get HF-vmoaos from dumpfile
  call rdsamo( ip2un, vmoao, norb, nbsf, rnel)
  !transform density from mo to ao basis
  Tvmoao = transpose( vmoao)
!  call tmtdag( pmo, nbsf, pao, nbsf, vmoao)
  call pmo2pao( pmo, nbsf, pao, nbsf, vmoao, .true.)

  open( 10, file = 'points', status = 'old', err = 222)
  rewind( 10)
  read(10,*) npnt

  write(6,'(/'' calulating density, vhartr and vcond at '',i6,'' points'')') npnt

  !loop over points
  do ipnt = 1, npnt

     vhartr = 0.d0
     vcond = 0.d0
!     vX = 0.d0
     rhovcond = 0.d0
!     rhovX = 0.d0

     read(10,*) Pcoord
     !calculate AO integrals at point Pcoord
     call cNATAOInt( Pcoord, NATAOInt, nbsf)
!     call cXAOInt(Pcoord, XAOInt, nbsf)
     !calculate the AO values at point Pcoord
     call cvalAO( Pcoord, valAO)
!     write(*,*) 'valAO'
!     write(*,*) valAO
!     write(*,*) 'pAO'
!     write(*,*) pAO
     !calculate the MO values at point Pcoord from AO values
     call vecmat( valAO, nbsf, Tvmoao, nbsf, valMO)
     !calculate hole density matrix
     call cond( rpcond, densty, norb, rnel, valMO, Pcoord, nbsf)
     !transform rpcond from MO to AO basis
!     call tmtdag( rpcond, nbsf, rpcndao, nbsf, vmoao)
     call pmo2pao( rpcond, nbsf, rpcndao, nbsf, vmoao, .true.)

     !calculate Hartree and conditional potential in AO basis
     do i = 1, nbsf

        ii = i * ( i - 1) / 2 

        do j = 1, i
           
           ij = ii + j
           vhartr = vhartr + pao(ij) * NATAOInt(ij)
!           vhartr = vhartr + pao(ij) * ( NATAOInt(i,j) + NATAOInt(j,i)) / 2
!           rhovX = rhovX + pao(ij) *( XAOInt(j,i) +  XAOInt(i,j))
           rhovcond = rhovcond + rpcndao(ij) * NATAOInt(ij)
!           rhovcond = rhovcond + rpcndao(ij) * ( NATAOInt(i,j) + NATAOInt(j,i)) / 2
        enddo
     enddo
     
!     vX = rhovX / densty
     vcond = rhovcond / densty

     write(11,'(3f10.3,e18.7)') Pcoord, vhartr / 2
     write(12,'(3f10.3,e18.7)') Pcoord, vcond
     write(13,'(3f10.3,e18.7)') Pcoord, densty
     write(14,'(3f10.3,e18.7)') Pcoord, rhovcond
     write(15,'(3f10.3,2e18.7)') Pcoord, ( vcond - vhartr) / 2, ( vcond - vhartr) / 2 * densty
!     write(16,'(3f10.3,f18.7)') Pcoord, vX

     call progress(ipnt, npnt)

  enddo

  write(6,'(''  ====>    done'')')

  close(20)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
!  close(16)
  write(*,*)

  return

222 write(6,'(''ERROR; file "points" does not exist'')')
  stop

end subroutine vbrain

subroutine vecmat(vec1, ndim1, xmat, ndim2, vec2)
  !***********************************************************************
  !
  ! multiplication vector * matrix
  !
  !
  !    vec2(i)  =  vec1(j) * xmat(j,i) sum over j
  !
  !
  !***********************************************************************
  implicit none

  integer, intent(in)           :: ndim1, ndim2
  double precision, intent(in)  :: vec1(ndim1), xmat(ndim2,ndim1)
  double precision, intent(out) :: vec2(ndim2)

  integer :: i, j

  vec2 = 0.d0

  do i = 1, ndim1
     do j = 1, ndim2

        vec2(j) = vec2(j) + vec1(i) * xmat(j,i)

     enddo
  enddo

  return
  
end subroutine vecmat

subroutine tmtdag( a, na, b, nb, q)
  !**********************************************************************
  !
  ! subroutine for matrix transformation : b = qaq+
  !
  !      b(kl) = q(ki)*a(ij)*q(lj)  sum over i,j
  !
  ! a  : symmetric matrix with dimension na*na stored in tridiagnal form
  !      and folded, i.e. a(i,j) = a_ij + a_ji
  !      element a(i,j) stored in position i*(i-1)/2+j (i>=j)
  ! b  : symmetric matrix with dimension nb*nb stored in tridiagnal form
  !      and folded, i.e. b(i,j) = b_ij + b_ji
  !      element b(i,j) stored in position i*(i-1)/2+j (i>=j)
  ! q  : transformation matrix with dimension nb*na
  !
  ! before transformation unfolding
  ! after transformation folding
  !
  !**********************************************************************
  implicit none

  integer, intent(in)           :: na, nb
  double precision, intent(in)  :: q(na,nb)
  double precision              :: a(na*(na+1)/2)
  double precision, intent(out) :: b(nb*(nb+1)/2)

  integer          :: i, ii, j, jj, k, l, kl
  double precision :: t(na)

  b = 0.d0
  
  !unfold a
  a = a * 0.5d0

  do i = 1, na

     ii = i * ( i + 1) / 2
     a(ii) = a(ii) * 2.d0

  enddo

  !multiplication
  do k = 1, nb
     do j = 1, na

        t(j) = 0.d0
        jj = j * ( j - 1) / 2

        do i = 1, j
           t(j) = t(j) + a(jj+i) * q(i,k)
        enddo

        do i = j + 1, na
           t(j) = t(j) + a(i*(i-1)/2+j) * q(i,k)
        enddo
     enddo

     do l = 1, k

        kl = k * ( k - 1) / 2 + l

        do j = 1, na
           b(kl) = b(kl) + t(j) * q(j,l)
        enddo
     enddo
  enddo

  !folding 
  b = b * 2.d0

  do i = 1, nb

     ii = i * ( i + 1) / 2
     b(ii) = b(ii) * 0.5d0

  enddo

  a = a * 2.d0

  do i = 1, na

     ii = i * ( i + 1) / 2
     a(ii) = a(ii) * 0.5d0

  enddo

  return

end subroutine tmtdag

subroutine pmo2pao( pmo, norb, pao, nbsf, vmoao, fold)
  !**********************************************************************
  !
  ! transforms density from MO to AO basis
  ! 
  !
  !   p_ao(k,l)  =  sum_mu,nu c_mu^k c_nu^l p_mo(mu,nu)
  !
  !
  ! with density stored in triangular form
  !
  !**********************************************************************
  implicit none

  integer, intent(in)           :: norb, nbsf
  double precision, intent(in)  :: pmo(norb*(norb+1)/2), vmoao(norb,nbsf)
  double precision, intent(out) :: pao(nbsf*(nbsf+1)/2)
  logical, intent(in)           :: fold

  integer          :: i, j, ii ,ij  
  double precision :: tpmo(norb,norb), tpao(nbsf,nbsf), tmat(nbsf,norb), tvmoao(norb,nbsf), ttpmo(norb*(norb+1)/2)
  double precision :: wrk(norb)
  
  tpao = 0.d0
  tpmo =  0.d0
  ttpmo = 0.d0

  if ( fold) then

     !transformation of pmo from triangluar stored to rectangular stored matrix
     !and unfolding
     do i = 1, norb

        ii = i * ( i - 1) / 2

        do j = 1, i

           ij = ii + j
           !tpmo(i,j) = pmo(ij)
           if( i .eq. j) then
              tpmo(i,i) = pmo(ij)
              ttpmo(ij) = pmo(ij)
           else
              tpmo(j,i) = pmo(ij) / 2
              tpmo(i,j) = tpmo(j,i)
              ttpmo(ij) = pmo(ij) / 2
           end if

        enddo

     enddo     

  else

     !transformation of pmo from triangluar stored to rectangular stored matrix
     do i = 1, norb

        ii = i * ( i - 1) / 2

        do j = 1, i

           ij = ii + j
           tpmo(j,i) = pmo(ij)
           tpmo(i,j) = pmo(ij)

        enddo

     enddo

  endif

  !matrix multiplication
!  call tftri( pao, ttpmo, tvmoao, wrk, nbsf, nbsf, nbsf)
  tvmoao = vmoao
!  tpao = matmul(matmul(tvmoao,tpmo),transpose(tvmoao))     
!  tpao = matmul(matmul(transpose(tvmoao),tpmo),tvmoao)
  call dsymm('r', 'u', nbsf, norb, 1.0d0, tpmo, norb, vmoao, nbsf, &
       & 0.0d0, tmat, nbsf)
!  call dsymm('l', 'u', nbsf, norb, 1.0d0, tpmo, norb, vmoao, nbsf, &
!       & 0.0d0, tmat, nbsf)
  call dgemm('n','t', nbsf, nbsf, norb, 1.0d0, tmat, nbsf, vmoao, &
       &  nbsf, 0.0d0, tpao, nbsf )

  !backtransform from rectangular to triangular saved matrix
  ! and folding
  if ( fold) then

     do i = 1, nbsf

        ii = i * ( i - 1) / 2

        do j = 1, i

           if ( abs( tpao(i,j) - tpao(j,i)) .gt. 1.d-8) then
              write(6,'(''ERROR IN MATRIX MULTIPLICTION'',2i4,2f22.6)') i, j, tpao(i,j), tpao(j,i)
              write(6,'(''subroutine pmo2pao'')')
              stop
           endif

           ij = ii + j

           if( i .eq. j) then
              pao(ij) = tpao(j,i)
           else
              pao(ij) = tpao(j,i) * 2
!              pao(ij) = pao(ij) * 2
           end if

        enddo

     enddo

  else

     do i = 1, nbsf

        ii = i * ( i - 1) / 2

        do j = 1, i

           if ( abs( tpao(i,j) - tpao(j,i)) .gt. 1.d-8) then
              write(6,'(''ERROR IN MATRIX MULTIPLICTION'',2i4,2f22.6)') i, j, tpao(i,j), tpao(j,i)
              write(6,'(''subroutine pmo2pao'')')
              stop
           endif

           ij = ii + j
           pao(ij) = tpao(j,i)

        enddo

     enddo

  endif

  return

end subroutine pmo2pao

subroutine tftri(h,f,t,wrk,m,n,ldt)
  !***************************************************************************
  !
  !     ----- transform the triangular matrix f using vectors t -----
  !                      h = t-dagger * f * t
  !     the order of the triangular matrix -h- is -m-   and   -f- is -n-
  !
  !                      (copy from gamess-us)
  !***************************************************************************
  implicit double precision(a-h,o-z)
  dimension h(*),f(*),t(ldt,m),wrk(n)
  parameter (zero=0.0d+00, one=1.0d+00, small=1.0d-11)

  m2 = ( m * m + m) / 2

  !        the computation here is h = t-dagger * (f * t),
  !        with the -dspmv- first producing one column of f*t,
  !        then the -dgemv- generates an entire row -j- of -h-.
  do j = 1, m

     ij = ( j * j - j) / 2
     call dspmv( 'u', n, one, f, t(1,j), 1, zero, wrk, 1)
     call dgemv( 't', n, j, one, t, ldt, wrk, 1, zero, h(ij+1), 1)

     do i=1,j
        if (abs(h(ij+i)).lt.small) h(ij+i)=zero
     enddo

  end do

  return

end subroutine tftri

subroutine progress(j, t)

  use ifport

  implicit none

  integer(kind=4), intent(in) :: j,t

  integer           :: k  
  character(len=37) :: bar="???% |                              |"  

  write(unit=bar(1:3),fmt="(i3)") 100 * j / t

  do k=1, int( dble(j) / t * 30)
     bar(6+k:6+k)="="
  enddo

  ! print the progress bar.                                                                                                                                                              
  write(*,fmt="(a1,a1,a37,$)") '+',char(13), bar

  return

end subroutine progress

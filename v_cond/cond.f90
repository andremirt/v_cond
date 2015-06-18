subroutine cond( rpcond, p1ex, norb, rnel, valMO, P, nbsf)
  !***************************************************************************
  !
  !  calculates hole density matrix from MO's
  !
  !
  !        rpcond(P)   =   sum_kl pdfmo_ijkl * f_k(P) * f_l(P)
  !
  !  
  !  the density p1(P) is calculated by
  !
  !
  !            p1(P)   =   sum_ij pmo_ij * f_i(P) * f_j(P)
  !
  !
  !  reference electron is in point P, the f's are the MO values at P
  !
  !  rpcond matrix in array 
  !  rpcond (in triangular form, rpcond(i*(i-1)/2+j) = sum of off-diagonal 
  !  elements (i,j) and (j,i) ). pcond matrix is calculated in same basis 
  !  as coef-matrix (usually mo's).
  !
  !  normalization :  pdfmo normalized to n(n-1)
  !                   pmo normalized to n
  !                   rpcond / p1 normalized to n-1
  !
  !  further arguments :  norb  : number of orbitals
  !                       rnel  : number of electrons
  !                       valmo : array to hold mo-values at P
  !
  !***************************************************************************
  use density
  implicit none

  integer, intent(in)           :: norb, nbsf
  double precision, intent(in)  :: rnel, valMO(nbsf), P(3)
  double precision, intent(out) :: rpcond(nbsf*(nbsf+1)/2), p1ex

  integer                    :: ntot, ntot2, m, iorb, i, j, ii, ij, k, l, kl
  double precision           :: tol, valij, valkl, cnel, tmp
  double precision,parameter :: explim = 18

  tol = dexp( -explim)
  ntot = norb * ( norb + 1) / 2
  ntot2 = ntot * ( ntot + 1) / 2

  p1ex = 0.d0

!!$  write(*,*) 'pmo'
!!$  write(*,*) pmo
!!$  write(*,*) 'valMO'
!!$  write(*,*) valMO

  do i = 1, nbsf

     ii = i * ( i - 1) / 2
     tmp = 0.d0

     do j = 1, i 

        ij = ii + j
        tmp = tmp + pmo(ij) * valMO(j)

     end do

     p1ex = p1ex + tmp * valMO(i)

  end do
  
  rpcond = 0.d0 

  do m = 1, ntot2
     
     i = ijkl(1,m)
     j = ijkl(2,m)
     k = ijkl(3,m)
     l = ijkl(4,m)
     ij = i * ( i - 1) / 2 + j
     kl = k * ( k - 1) / 2 + l
     valij = pdfmo(m) * valMO(i) * valMO(j)
     valkl = pdfmo(m) * valMO(k) * valMO(l)
     rpcond(ij) = rpcond(ij) + valkl
     rpcond(kl) = rpcond(kl) + valij

  enddo

  if ( p1ex .lt. tol) then
     write(6,'(''WARNING: small density in reference position'')')
     write(6,'(''  rho : '',g16.8,''  point : '',3g11.3)') p1ex, P
  endif

  cnel=0.d0

  do iorb = 1 , norb
     cnel = cnel + rpcond(iorb*(iorb+1)/2)
  enddo

  cnel = cnel / p1ex 

  if ( abs( cnel + 1.d0 - rnel) .gt. 1.d-6) then
     write(6,'(''ERROR; cnel = '',f16.10,'' rnel = '',f16.10, &
          & '' diff = '',g16.8)') cnel, rnel, cnel + 1.d0 - rnel
     stop
  endif

  return

end subroutine cond

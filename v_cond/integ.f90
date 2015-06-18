module integral
  !********************************************************************
  !
  !  module containing routines in relation 
  !    to nuclear attraction type integral over AO's
  !
  !                    
  !                        /
  !      I(k,l,r1)   =  -  | f_k(Ra-r2)f_l(Rb-r2) / |r1-r2| dr2
  !                        /
  !                   
  !
  !  with f_k, f_l contracted normalized cartesian gaussians
  !    cenetered at Ra resp. Rb and nuclear charge replaced by -1
  !
  !  - spherical part of gaussians shall be given cartesian
  !
  !********************************************************************
  use readgauss
  implicit none

  private
  public cNATAOInt, cOlAO, boysf
  
  integer                     :: allst
  double precision, external  :: dfac, dbinom

  type tboysf
     integer          :: nu
     double precision :: u
  end type tboysf
  
  type(tboysf) :: boysf

contains

  subroutine cNATAOInt( point, NATAOInt, nbsf)
    !******************************************************************
    !
    ! main subroutine for AO-integral calculation in point r1
    !
    ! - contains loop over all AO combinations
    ! - loop over primitive gaussians
    !
    !******************************************************************
    integer, intent(in)           :: nbsf
    double precision, intent(in)  :: point(3)
    double precision, intent(out) :: NATAOInt(nbsf*(nbsf+1)/2)
!    double precision, intent(out) :: NATAOInt(nbsf,nbsf)

    integer                       :: at1, at2, sh1, sh2, pr1, pr2, LA, LB, MA, MB, &
         & lmn1(3), lmn2(3), isp1, isp2, iorb1, iorb2, indx, ubsh2, ubMB, i1i1
    double precision              :: alpha, beta, centr1(3), centr2(3), coeff1, coeff2
    
    NATAOInt = 0.0d0
    iorb1 = 0
    
    do at1 = 1, nat

       centr1 = geom(at1)%coord

       do sh1 = 1, bset(at1)%nshells

          LA = bset(at1)%primtvs(sh1)%l
!!$          if ( LA .eq. 8) LA = 5

          do MA = 1, AngMTab(LA)%ncomb

             lmn1(:) = AngMTab(LA)%lmn(:,MA)
             iorb1 = iorb1 + 1
             i1i1 = iorb1 * ( iorb1 - 1) / 2 
             iorb2 = 0

             do at2 = 1, at1
!             do at2 = 1, nat

                centr2 = geom(at2)%coord

                if ( at1 .eq. at2) then
                   ubsh2 = sh1
                else
                   ubsh2 = bset(at2)%nshells
                end if
                
                do sh2 = 1, ubsh2

                   LB = bset(at2)%primtvs(sh2)%l
!!$                   if ( LB .eq. 8) LB = 5

                   if ( at1 .eq. at2 .and. sh1 .eq. sh2) then
                      ubMB = MA
                   else
                      ubMB = AngMTab(LB)%ncomb
                   end if

                   do MB = 1, ubMB

                      lmn2(:) = AngMTab(LB)%lmn(:,MB)
                      iorb2 = iorb2 + 1
                      indx = i1i1 + iorb2

                      do pr1 = 1, bset(at1)%primtvs(sh1)%nprims

                         coeff1 = bset(at1)%primtvs(sh1)%NormCoeff(pr1,MA)
                         alpha = bset(at1)%primtvs(sh1)%exp(pr1)

                         do pr2 = 1, bset(at2)%primtvs(sh2)%nprims
                
                            coeff2 = bset(at2)%primtvs(sh2)%NormCoeff(pr2,MB)
                            beta = bset(at2)%primtvs(sh2)%exp(pr2)

!!$                            write(*,*) 'here2', at1, at2, iorb1, iorb2
!!$                            write(*,*) 'sh1', sh1, MA, sh2, MB
                            NATAOInt(indx) = NATAOInt(indx) + coeff1 * coeff2 * &
                                  & NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!                            NATAOInt(iorb1,iorb2) = NATAOInt(iorb1,iorb2) + coeff1 * coeff2 * &
!                                  & NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!                            coeff2 = NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!                            write(*,*) 'NATInt', iorb1, iorb2, pr1, pr2, coeff2
!!$
!!$                            if ( LA .eq. 5) then !sp-shell on at1
!!$                               
!!$                               coeff1 = bset(at1)%primtvs(sh1)%coeff2(pr1)
!!$
!!$                               do isp1 = 1, 3
!!$                
!!$                                  lmn1 = AngMTab(LA)%lmn(:,isp1+1)
!!$                                  indx = ( iorb1 + isp1) * ( iorb1 + isp1 - 1) / 2 + iorb2
!!$                                  NATAOInt(indx) = NATAOInt(indx) + coeff1 * coeff2 * &
!!$                                       & NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!!$
!!$                                  if ( LB .eq. 5) then !sp-shell on at1 and at2
!!$
!!$                                     coeff2 = bset(at2)%primtvs(sh2)%coeff2(pr2)
!!$
!!$                                     do isp2 = 1, 3
!!$                
!!$                                        lmn2 = AngMTab(LB)%lmn(:,isp2+1)
!!$                                        indx = ( iorb1 + isp1) * ( iorb1 + isp1 - 1) / 2 + iorb2 + isp2
!!$                                        NATAOInt(indx) = NATAOInt(indx) + coeff1 * coeff2 * &
!!$                                             & NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!!$
!!$                                     end do
!!$                               
!!$                                     iorb2 = iorb2 + 3
!!$
!!$                                  end if
!!$                                  
!!$                               end do
!!$                   
!!$                               iorb1 = iorb1 + 3
!!$                               
!!$                            end if
!!$                            
!!$                            if ( LB .eq. 5) then !sp-shell on at2
!!$
!!$                               coeff2 = bset(at2)%primtvs(sh2)%coeff2(pr2)
!!$
!!$                               do isp2 = 1, 3
!!$
!!$                                  lmn2 = AngMTab(LB)%lmn(:,isp2+1)
!!$                                  indx = iorb1 * ( iorb1 - 1) / 2 + iorb2 + isp2
!!$                                  NATAOInt(indx) = NATAOInt(indx) + coeff1 * coeff2 * &
!!$                                       & NATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)
!!$
!!$                               end do
!!$
!!$                               iorb2 = iorb2 + 3
!!$
!!$                            end if

                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    
    return

  end subroutine cNATAOInt

  subroutine cOlAO
    !********************************************************************************
    !
    ! calculates the onside overlap integral over contracted gaussians 
    !  and normalizes the AOs resp. gaussians
    !
    !********************************************************************************
    integer          :: iat, ish, il, ipr1, ipr2, lmn(3)
    double precision :: OlInt, alpha, beta, coeff1, coeff2

    write(6,'('' calculating onside overlap for all AOs to check normalization'')')

    do iat = 1, nat

       do ish = 1, bset(iat)%nshells

          do il = 1, bset(iat)%primtvs(ish)%l

             OlInt = 0.d0
             lmn = AngMTab(bset(iat)%primtvs(ish)%l)%lmn(:,il)

             do ipr1 = 1, bset(iat)%primtvs(ish)%nprims

                alpha = bset(iat)%primtvs(ish)%exp(ipr1)
                coeff1 = bset(iat)%primtvs(ish)%NormCoeff(ipr1,il)

                do ipr2 = 1, bset(iat)%primtvs(ish)%nprims

                   beta = bset(iat)%primtvs(ish)%exp(ipr2)
                   coeff2 = bset(iat)%primtvs(ish)%NormCoeff(ipr2,il)
                   OlInt = OlInt + coeff1 * coeff2 * OlPInt( lmn, alpha, beta)

                end do
             end do

             if( abs(1.d0 - OlInt) .gt. 1.d-8) then
                write(6,*) 'ERROR: something went wrong with normalization'
                write(6,*) 'diagonalelement of overlap matrix not equal 1'
                write(6,'(''Atom : '',i4,''   shell : '',i4)') iat, ish
                write(6,'(''S(AO_'',i1,'','',i1,''|AO_'',i1,'','',i1,'') =  '',f20.6)')  iat, ish, iat, ish, OlInt
                stop
             end if

          end do
       end do
    end do

    write(6,*) '   ====>    all AOs are proper normalized to 1'

    return

  end subroutine cOlAO

  double precision function NATInt( lmn1, lmn2, alpha, beta, A, B, R1)
    !*********************************************************************************
    !
    ! calculates the nuclear attraction type integral
    !
    ! 
    !                      /
    !     NATInt(R2)   =   | d3R1 pgA_lmn1(R1) * pgB_lmn2(R1) / |R1 - R2|
    !                      /
    !                      
    !
    !   with nuclear charge replaced by -1
    !
    ! pg_lmnA is a primitive cartesian gaussian centered at A
    !
    !
    !     pgA_lmn1   =   (x-Ax)**l1 * (y-Ay)**m1 * (z-Az)**n1 * exp( -a(R2-A)**2)
    !
    !
    ! this routine contains: 
    !   1st step: Gaussian product theorem is applied
    !   2nd step: prefactors of form  (x - A)**l are expressed by 
    !               Hermite polynomials
    !
    ! e.g.: Petersson, Hellsing, EurJPhys 31 (2010), p.37
    !
    !********************************************************************************
    integer, intent(in)          :: lmn1(3), lmn2(3)
    double precision, intent(in) :: alpha, beta, A(3), B(3), R1(3)

    integer          :: coord, l1, l2, m1, m2, n1, n2, lfac1, lfac2, &
         & mfac1, mfac2, nfac1, nfac2
    double precision :: g, n, P(3), mnsAmBsq, AmB(3), PmR1(3), mgsPmR1sq, &
         & sum1, sum2, sum3, sum4, sum5, sum6

!!$    write(*,*) 'NATInt, lmn1, lmn2, alpha, beta, A, B, R1'
!!$    write(*,*) lmn1
!!$    write(*,*) lmn2
!!$    write(*,*) alpha, beta
!!$    write(*,*) A
!!$    write(*,*) B
!!$    write(*,*) R1

    g = alpha + beta
    n = alpha * beta / g
    P = ( alpha * A + beta * B) / g
    AmB = A - B
    mnsAmBsq = - n * sum( AmB * AmB)
    PmR1 = P - R1
    mgsPmR1sq = - g * sum( PmR1 * PmR1)
    sum1 = 0.0d0

    do l1 = 0, int( lmn1(1) / 2.)

       lfac1 = lmn1(1) - 2 * l1
       sum2 = 0.d0

       do l2 = 0, int( lmn2(1) / 2.)

          lfac2 = lmn2(1) - 2 * l2
          sum3 = 0.d0

          do m1 = 0, int( lmn1(2) / 2.)

             mfac1 = lmn1(2) - 2 * m1
             sum4 = 0.d0
             
             do m2 = 0, int( lmn2(2) / 2.)

                mfac2 = lmn2(2) - 2 * m2
                sum5 = 0.d0

                do n1 = 0, int( lmn1(3) / 2.)

                   nfac1 = lmn1(3) - 2 * n1
                   sum6 = 0.d0

                   do n2 = 0, int( lmn2(3) / 2.)

                      nfac2 = lmn2(3) - 2 * n2             
                      sum6 = sum6 + derivExpAll( lfac1, lfac2, mfac1, mfac2, nfac1, nfac2, n, AmB, g, &
                           & PmR1, mgsPmR1sq, alpha, beta) / ( dfac( n2) * dfac( nfac2) &
                           & * beta**( lmn2(3) - n2))

                   end do

                   sum5 = sum5 + sum6 / ( dfac( n1) * dfac( nfac1) &
                        & * alpha**( lmn1(3) - n1))

                end do
                
                sum4 = sum4 + sum5 / ( dfac( m2) * dfac( mfac2) &
                     & * beta**( lmn2(2) - m2))

             end do

             sum3 = sum3 + sum4 / ( dfac( m1) * dfac( mfac1) &
                  & * alpha**( lmn1(2) - m1))

          end do
             
          sum2 = sum2 + sum3 / ( dfac( l2) * dfac( lfac2) &
               & * beta**( lmn2(1) - l2))

       end do

       sum1 = sum1 + sum2 / ( dfac( l1) * dfac( lfac1) &
            & * alpha**( lmn1(1) - l1))
       
    end do

    NATInt = sum1 * dfac( lmn1(1)) * dfac( lmn2(1)) * dfac( lmn1(2)) * &
         & dfac( lmn2(2)) * dfac( lmn1(3)) * dfac( lmn2(3)) / &
         & 2**( sum( lmn1) + sum( lmn2)) * 2 * Pi / g * dexp( mnsAmBsq)

!!$    write(*,*) 'result', NATInt

    return

  end function NATInt
  
  double precision function OlPInt( lmn, alpha, beta)
    !******************************************************************************
    !
    ! calculates overlap integral of primitive gaussian
    !
    ! 
    !                  /
    !     OlPInt   =   | d3R1 pg_lmn(R1) * pg_lmn(R1)
    !                  /
    !                      
    ! this routine contains: 
    !   prefactors of form x**l are expressed by Hermite polynomials
    !
    ! e.g.: Petersson, Hellsing, EurJPhys 31 (2010), p.37
    !
    !******************************************************************************
    integer, intent(in)          :: lmn(3)
    double precision, intent(in) :: alpha, beta

    integer          :: iL, iM, iN, il1, il2, im1, im2, in1, in2, &
         & ltwoil1, mtwoim1, ntwoin1, ltwoil2, mtwoim2, ntwoin2
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, g, nu

    g = alpha + beta
    nu = alpha * beta / g
    sum1 = 0.0d0

    do il1 = 0, int( lmn(1) / 2.)

       sum2 = 0.d0
       ltwoil1 = lmn(1) - 2 * il1

       do il2 = 0, int( lmn(1) / 2.)

          sum3 = 0.d0
          ltwoil2 = lmn(1) - 2 * il2
          
          do im1 = 0, int( lmn(2) / 2.)
       
             sum4 = 0.d0
             mtwoim1 = lmn(2) - 2 * im1
          
             do im2 = 0, int( lmn(2) / 2.)
               
                sum5 = 0.d0
                mtwoim2 = lmn(2) - 2 * im2
          
                do in1 = 0, int( lmn(3) / 2.)

                   sum6 = 0.d0
                   ntwoin1 = lmn(3) - 2 * in1
          
                   do in2 = 0, int( lmn(3) / 2.)

                      ntwoin2 = lmn(3) - 2 * in2
                      iN = ntwoin1 + ntwoin2
                      sum6 = sum6 + dfac( iN) * derivExpAB( iN, nu, 0.d0) &
                           & / ( dfac( in2) * dfac( ntwoin2) * beta**( lmn(3) - in2))

                   end do
                   
                   sum5 = sum5 + sum6 / ( dfac( in1) * dfac( ntwoin1) &
                           & * alpha**( lmn(3) - in1)) * ( -1)**(ntwoin1)

                end do
                
                iM = mtwoim1 + mtwoim2
                sum4 = sum4 + sum5 * dfac( iM) * derivExpAB( iM, nu, 0.d0) &
                     & / ( dfac( im2) * dfac( mtwoim2) * beta**( lmn(2) - im2))

             end do

             sum3 = sum3 + sum4 / ( dfac( im1) * dfac( mtwoim1) &
                  & * alpha**( lmn(2) - im1)) * ( -1)**(mtwoim1)
             
          end do

          iL = ltwoil1 + ltwoil2
          sum2 = sum2 + sum3 * dfac( iL) * derivExpAB( iL, nu, 0.d0) &
               & / ( dfac( il2) * dfac( ltwoil2) * beta**( lmn(1) - il2))

       end do

       sum1 = sum1 + sum2 / ( dfac( il1) * dfac( ltwoil1) &
            & * alpha**( lmn(1) - il1)) * ( -1)**(ltwoil1)

    end do
    
    OlPInt = sum1 * dfac( lmn(1)) * dfac( lmn(1)) * dfac( lmn(2)) * &
         & dfac( lmn(2)) * dfac( lmn(3)) * dfac( lmn(3)) / &
         & 2**( 2 * sum( lmn)) * dsqrt(( Pi / g)**3)

    return

  end function  OlPInt

  double precision function derivExpAll( l1, l2, m1, m2, n1, n2, nu, AmB, g, PmR1, mgsPmR1sq, a, b)
    !******************************************************************************
    !
    ! 3rd step: Leibnitz theorem is applied for the six derivatives
    !
    !                                   
    !                    d**l1   d**l2    d**m1   d**m2    d**n1   d**n2   
    !   derivExpPR1  =  ------- -------  ------- -------  ------- ------- *
    !                   dAx**l1 dBx**l2  dAy**m1 dBy**m2  dAz**n1 dBz**n2 
    !                                   
    !                                           1
    !                                           /
    !                   * exp( -n * (A-B)**2) * | dt exp( -g(P - R1)**2 * t**2)
    !                                           /
    !                                           0
    !
    !******************************************************************************
    integer, intent(in)          :: l1, l2, m1, m2, n1, n2
    double precision, intent(in) :: nu, AmB(3), g, PmR1(3), mgsPmR1sq, a, b

    integer          :: L, M, N, iL, iM, iN, il1, il2, im1, im2, in1, in2
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6

!    L = l1 + l2
!    M = m1 + m2
!    N = n1 + n2

    sum1 = 0.0d0

    do il1 = 0, l1

       sum2 = 0.d0

       do il2 = 0, l2

          sum3 = 0.d0
          
          do im1 = 0, m1
       
             sum4 = 0.d0
          
             do im2 = 0, m2
               
                sum5 = 0.d0
          
                do in1 = 0, n1

                   sum6 = 0.d0
          
                   do in2 = 0, n2

                      iN = in1 + in2
                      sum6 = sum6 + dbinom( n2, in2) * derivExpAB( iN, nu, AmB(3)) * dfac( iN) * &
                           & derivExpPR1( l1 - il1, l2 - il2, m1 - im1, m2 - im2, n1 - in1, n2 - in2, &
                           & g, PmR1, mgsPmR1sq, a, b)

                   end do
       
                   sum5 = sum5 + dbinom( n1, in1) * sum6 * ( -1)**(in1)

                end do
       
                iM = im1 + im2
                sum4 = sum4 + dbinom( m2, im2) * derivExpAB( iM, nu, AmB(2)) * sum5 * dfac( iM)

             end do
    
             sum3 = sum3 + dbinom( m1, im1) * sum4 * ( -1)**(im1)

          end do
       
          iL = il1 + il2
          sum2 = sum2 + dbinom( l2, il2) * derivExpAB( iL, nu, AmB(1)) * sum3 * dfac( iL)

       end do

       sum1 = sum1 + dbinom( l1, il1) * sum2 * ( -1)**(il1)

    end do

    derivExpAll = sum1

    return

  end function derivExpAll
  
  double precision function derivExpAB( j, n, AmB)
    !***********************************************************************
    !
    ! 4th step: calculates the derivative
    !
    !                                   
    !                    d**j1    d**j2   
    !   derivExpPR1  =  -------  ------- exp( -n (A - B)**2)
    !                    dAc**j1 dBc**j2 
    !                                   
    !
    !   with c = {x,y,z} and j = j1 + j2
    !
    ! the factors ( -1)**j1 and j! are in routine derivExpAll
    ! the factor exp( -n (A - B)**2) is in routine NATInt
    !
    !***********************************************************************
    integer, intent(in)          :: j
    double precision, intent(in) :: n, AmB

    integer          :: r, jmtwor
    double precision :: sum

    sum = 0.0d0

    do r = 0, int( j / 2.)

       jmtwor = j - 2 * r 
       sum = sum + ( -1)**r * 2**( jmtwor) * n**( j - r) * AmB**( jmtwor) / &
            & ( dfac( r) * dfac( jmtwor))

    end do
    
    derivExpAB = sum

    return

  end function derivExpAB

  double precision function derivExpPR1( l1, l2, m1, m2, n1, n2, g, PmR1, mgsPmR1sq, a, b)
    !*******************************************************************************
    !
    ! 5th step: calculates the derivative of the integral
    !
    !
    !                    d**l1   d**l2    d**m1   d**m2    d**n1   d**n2   
    !   derivExpPR1  =  ------- -------  ------- -------  ------- ------- *
    !                   dAx**l1 dBx**l2  dAy**m1 dBy**m2  dAz**n1 dBz**n2 
    !
    !
    !                     1
    !                     /
    !                   * | dt exp( -g(P - R1)**2 * t**2)
    !                     /
    !                     0
    !
    !
    !   by means of boys integral
    !
    !*******************************************************************************
    integer, intent(in)          :: l1, l2, m1, m2, n1, n2
    double precision, intent(in) :: a, b, g, PmR1(3), mgsPmR1sq

    integer          :: s, t, u, L, M, N, Lmtwos, Mmtwot, Nmtwou
    double precision :: sum1, sum2, sum3

    L = l1 + l2
    M = m1 + m2
    N = n1 + n2
    sum1 = 0.0d0

    do s = 0, int( L / 2.)

       Lmtwos =  L - 2 * s
       sum2 = 0.0d0

       do t = 0, int( M / 2.)

          Mmtwot =  M - 2 * t
          sum3 = 0.0d0

          do u = 0, int( N / 2.)

             Nmtwou =  N - 2 * u
             sum3 = sum3 + ( -1)**( N + u) * 2**(Nmtwou) * PmR1(3)**( Nmtwou) / &
                  & ( dfac( u) * dfac( Nmtwou) * g**u) * Fnu( mgsPmR1sq, L + M + N - s - t - u)

          end do

          sum2 = sum2 + ( -1)**( M + t) * 2**(Mmtwot) * PmR1(2)**( Mmtwot) / &
               & ( dfac( t) * dfac( Mmtwot) * g**t) * sum3

       end do

       sum1 = sum1 + ( -1)**( L + s) * 2**(Lmtwos) * PmR1(1)**( Lmtwos) / &
            & ( dfac( s) * dfac( Lmtwos) * g**s) * sum2

    end do
    
    derivExpPR1 = a**( l1 + m1 + n1) * b**( l2 + m2 + n2) * dfac( L) * dfac( M) * dfac( N) * sum1

    return

  end function derivExpPR1

  recursive double precision function Fnu( u, nu)
    !***********************************************************************
    !
    ! 6th step: calculates the boys integral
    !
    !                  1
    !                  /
    !    F_nu(u)   =   | dt t**( 2 * u) * exp( -u * t**2)
    !                  /
    !                  0
    !
    !  numericaly by help of dqag routine of SLATEC
    !
    !***********************************************************************
    integer, intent(in)          :: nu
    double precision, intent(in) :: u

    integer                    :: k, neval, ier, last, iwork(100)
    double precision           :: sum, vresult, verror, work(1000)
    double precision, external :: boys

    !possibility 1: numerical integration
    boysf%nu = 2 * nu
    boysf%u = u

    call dqag( boys, 0.d0, 1.d0, 1d-14, 1.0d-11, 1, vresult, verror, neval, ier, &
         & 100, 1000, last, iwork, work)

    Fnu = vresult

    return

!!$    this ways are not working if position of ``nuclear charge'' 
!!$    close to origin of gaussian
!!$
!!$    !possibility 2
!!$    Fnu = dgami(nu + 0.5d0, u) / ( 2 * u**(nu + 0.5))
!!$    write(*,*) 'Fnu0', u, nu, Fnu
!!$    return
!!$
!!$    !possibility 3
!!$    if ( nu .eq. 0) then
!!$       !use downward recursion
!!$       Fnu = ( 2 * u * Fnu( u, nu + 1) + dexp( -u)) / ( 2 * nu + 1)
!!$       write(*,*) 'Fnu1', u, nu, Fnu
!!$       return
!!$
!!$    end if
!!$    
!!$    if( abs( u**nu) .lt. 1.d-28) then
!!$       !simplified integration with exp(- u * t**2) = 1
!!$       Fnu = 1.d0 / ( 2 * nu + 1)
!!$       write(*,*) 'tiny', tiny(1.0d0), tiny(1.0d0)**(1/nu)
!!$       write(*,*) 'Fnu2', u, nu, Fnu
!!$       return
!!$
!!$    end if
!!$
!!$    sum = 0.0d0
!!$    
!!$    do k = 0, nu - 1
!!$ 
!!$       sum = sum + dfac( nu - k) / ( 4**k * dfac( 2 * nu - 2 * k) * u**( k + 1))
!!$       write(*,*) 'sumF-nuloop', u, nu, k
!!$       write(*,*) 'sumF-nuloop', sum
!!$    
!!$    end do
!!$    
!!$    write(*,*) 'sumF-nu', u, nu, sum
!!$
!!$    Fnu = dfac( 2 * nu) / ( 2 * dfac( nu)) * ( dsqrt( Pi) / &
!!$         & ( 4**nu * u**( nu + 0.5d0)) * derf( dsqrt( u)) - sum * dexp( - u))
!!$
!!$    write(*,*) 'Fnu3', u, nu, Fnu
!!$
!!$    return

  end function Fnu

end module integral

double precision function boys( x)
  !****************************************************************
  !
  ! the boys function is defined as it is used in dqag routine 
  !  for numerical integration
  !
  !****************************************************************
  use integral
  implicit none
  double precision, intent(in) :: x

  boys = x**boysf%nu * exp( boysf%u * x**2)

  return

end function boys

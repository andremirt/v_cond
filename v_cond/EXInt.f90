!! unfinished subrouines on calculation of exchange part

  subroutine cXAOInt( point, XAOInt, nbsf)
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
    double precision, intent(out) :: XAOInt(nbsf,nbsf)

    integer                       :: at1, at2, sh1, sh2, pr1, pr2, LA, LB, MA, MB, &
         & lmn1(3), lmn2(3), isp1, isp2, iorb1, iorb2, indx, ubsh2, ubMB
    double precision              :: alpha, beta, centr1(3), centr2(3), coeff1, coeff2
    
    XAOInt = 0.0d0
    iorb1 = 0
    
    do at1 = 1, nat

       centr1 = geom(at1)%coord

       do sh1 = 1, bset(at1)%nshells

          LA = bset(at1)%primtvs(sh1)%l
!!$          if ( LA .eq. 8) LA = 5

          do MA = 1, AngMTab(LA)%ncomb

             lmn1(:) = AngMTab(LA)%lmn(:,MA)
             iorb1 = iorb1 + 1
             iorb2 = 0

             do at2 = 1, nat

                centr2 = geom(at2)%coord

!                if ( at1 .eq. at2) then
!                   ubsh2 = sh1
!                else
                   ubsh2 = bset(at2)%nshells
!                end if
                
                do sh2 = 1, ubsh2

                   LB = bset(at2)%primtvs(sh2)%l
!!$                   if ( LB .eq. 8) LB = 5

!                   if ( at1 .eq. at2 .and. sh1 .eq. sh2) then
!                      ubMB = MA
!                   else
                      ubMB = AngMTab(LB)%ncomb
!                   end if

                   do MB = 1, ubMB

                      lmn2(:) = AngMTab(LB)%lmn(:,MB)
                      iorb2 = iorb2 + 1
                      indx = iorb1 * ( iorb1 - 1) / 2 + iorb2

                      do pr1 = 1, bset(at1)%primtvs(sh1)%nprims

                         coeff1 = bset(at1)%primtvs(sh1)%NormCoeff(pr1,MA)
                         alpha = bset(at1)%primtvs(sh1)%exp(pr1)

                         do pr2 = 1, bset(at2)%primtvs(sh2)%nprims
                
                            coeff2 = bset(at2)%primtvs(sh2)%NormCoeff(pr2,MB)
                            beta = bset(at2)%primtvs(sh2)%exp(pr2)

!!$                            write(*,*) 'here2', at1, at2, iorb1, iorb2
!!$                            write(*,*) 'sh1', sh1, MA, sh2, MB
                            XAOInt(iorb1,iorb2) = XAOInt(iorb1,iorb2) + coeff1 * coeff2 * &
                                  & XATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)

                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
    
    return

  end subroutine cXAOInt

  double precision function XATInt( lmn1, lmn2, alpha, beta, centr1, centr2, point)

    integer, intent(in)          :: lmn1(3), lmn2(3)
    double precision, intent(in) :: alpha, beta, centr1(3), centr2(3), point(3)

    integer          :: l2, m2, n2, lfac2, mfac2, nfac2
    double precision :: prefac, sum1, sum2, sum3, CmR1(3), mbsCmR1sq


    prefac = 2 * ( point(1) - centr1(1))**lmn1(1) *( point(2) - centr1(2))**lmn1(2) * ( point(3) - centr1(3))**lmn1(3) * dexp( - alpha * sum( point - centr1) * sum( point - centr1)) / dsqrt( Pi)
    CmR1 = centr2 - point
    mbsCmR1sq = -beta * sum( CmR1 * CmR1)
    sum1= 0.d0

    do l2 = 0, int( lmn2(1) / 2.)

       lfac2 = lmn2(1) - 2 * l2
       sum2 = 0.d0

       do m2 = 0, int( lmn2(2) / 2.)

          mfac2 = lmn2(2) - 2 * m2
          sum3 = 0.d0

          do n2 = 0, int( lmn2(3) / 2.)

             nfac2 = lmn2(3) - 2 * n2
             sum3 = sum3 + derivExpPR1X( lfac2, mfac2, nfac2, &
                           & CmR1, beta, mbsCmR1sq) / ( dfac( n2) * dfac( nfac2) &
                           & * beta**( lmn2(3) - n2))

          end do
          
          sum2 = sum2 + sum3 / ( dfac( m2) * dfac( mfac2) &
               & * beta**( lmn2(2) - m2))

       end do
       
       sum1 = sum1 + sum2 / ( dfac( l2) * dfac( lfac2) &
            & * beta**( lmn2(1) - l2))
       
    end do
             
    XATInt = sum1 * dfac( lmn2(1)) * dfac( lmn2(2)) * &
         & dfac( lmn1(3)) * dfac( lmn2(3)) / &
         & 2**( sum( lmn2)) * 2 * Pi / beta

    return

  end function XATInt

!!$  double precision function derivExpAllX( l2, m2, n2, centr2, point, beta)
!!$
!!$    integer, intent(in)          :: l2, m2, n2
!!$    double precision, intent(in) :: centr(3), point(3), beta
!!$
!!$    integer          :: iL, iM, iN, il2, im2, in2
!!$    double precision :: sum1, sum2, sum3, AmB(3)
!!$
!!$    AmB = centr2 - point
!!$    sum1 = 0.0d0
!!$
!!$    do il2 = 0, l2
!!$
!!$       sum2 = 0.d0
!!$
!!$       do im2 = 0, m2
!!$
!!$          sum3 = 0.d0
!!$
!!$          do in2 = 0, n2
!!$
!!$             sum3 = sum3 + dbinom( n2, in2) * derivExpAB( in2, beta, AmB(3)) * dfac( in2) * &
!!$                  & derivExpPR1( l2 - il2, m2 - im2, n2 - in2,)
!!$
!!$          end do
!!$
!!$          sum2 = sum2 + dbinom( m2, im2) * derivExpAB( im2, beta, AmB(2)) * sum3 * dfac( im2)
!!$
!!$       end do
!!$
!!$       sum1 = sum1 + dbinom( l2, il2) * derivExpAB( il2, beta, AmB(3)) * sum2  * dfac( il2)
!!$
!!$    end do
!!$
!!$    derivExpAllX = sum1
!!$
!!$    return
!!$
!!$  end function derivExpAllX

  double precision function derivExpPR1X( l2, m2, n2, CmR1, beta, mbsCmR1sq)

    integer, intent(in)          :: l2, m2, n2
    double precision, intent(in) :: beta, CmR1(3), mbsCmR1sq

    integer          :: s, t, u, l2mtwos, m2mtwot, n2mtwou
    double precision :: sum1, sum2, sum3

    sum1 = 0.0d0

    do s = 0, int( l2 / 2.)

       l2mtwos =  l2 - 2 * s
       sum2 = 0.0d0

       do t = 0, int( m2 / 2.)

          m2mtwot =  m2 - 2 * t
          sum3 = 0.0d0

          do u = 0, int( n2 / 2.)

             n2mtwou =  n2 - 2 * u
             sum3 = sum3 + 2**(n2mtwou) * CmR1(3)**( n2mtwou) / &
                  & ( dfac( u) * dfac( n2mtwou) * beta**u) * Fnu( mbsCmR1sq, l2 + m2 + n2 - s - t - u) * ( -1)**( n2 + u)

          end do

          sum2 = sum2 + ( -1)**( m2 + t) * 2**(m2mtwot) * CmR1(2)**( m2mtwot) / &
               & ( dfac( t) * dfac( m2mtwot) * beta**t) * sum3

       end do

       sum1 = sum1 + ( -1)**( l2 + s) * 2**( l2mtwos) * CmR1(1)**( l2mtwos) / &
            & ( dfac( s) * dfac( l2mtwos) * beta**s) * sum2

    end do
    
    derivExpPR1X = beta**( l2 + m2 + n2) * dfac( l2) * dfac( m2) * dfac( n2) * sum1

    return

  end function derivExpPR1X

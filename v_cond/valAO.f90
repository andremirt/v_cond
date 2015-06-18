module mvalAO

  use readgauss
  implicit none
  public cvalAO

contains

  subroutine cvalAO( P, valAO)
    !************************************************************************
    !
    ! calculates the value of the AO's at point P, where the AO's are 
    !  normalized contracted gaussians
    !
    !************************************************************************
    double precision, intent(in)  :: P(3)
    double precision, intent(out) :: valAO(:)

    integer          :: iat, ish, ipr, iam, l, indx
    double precision :: center(3), cexpp, prefacp

    valAO = 0.d0
    indx = 0
!    write(*,*) 'valAO enter'
    do iat =  1, nat

       center = geom(iat)%coord

       do ish = 1, bset(iat)%nshells

          l = bset(iat)%primtvs(ish)%l
 
          do iam = 1, AngMTab(l)%ncomb

             cexpp = 0.0d0

             do ipr = 1, bset(iat)%primtvs(ish)%nprims

                cexpp = cexpp + bset(iat)%primtvs(ish)%NormCoeff(ipr,iam) * &
                     & dexp( - bset(iat)%primtvs(ish)%exp(ipr) * &
!                     & sum( (center - P) * (center - P)))
                     & ((P(1)-center(1))**2 + (P(2)-center(2))**2 + (P(3)-center(3))**2))
             end do
!    write(*,*) 'valAO mid'

!    write(*,*) 'valAO mid'
             indx = indx + 1
             prefacp = ( P(1) - center(1))**AngMTab(l)%lmn(1,iam) * &
                  & ( P(2) - center(2))**AngMTab(l)%lmn(2,iam) * &
                  & ( P(3) - center(3))**AngMTab(l)%lmn(3,iam)

             valAO(indx) = prefacp * cexpp

!             write(*,*) 'valAO', indx, iat, ish, l, iam, valAO(indx)
          
          end do

       end do
    end do

    return

  end subroutine cvalAO

end module mvalAO

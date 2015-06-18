subroutine rdsamo( idmp, vmoao, norb, nbsf, rnel)
  !************************************************************************
  !
  ! Read HF-MOs in AO basis from the dumpfile and store in vmoao
  !
  !************************************************************************
  implicit none

  integer, intent(in)           :: idmp, norb, nbsf
  double precision, intent(in)  :: rnel
  double precision, intent(out) :: vmoao(nbsf,nbsf)

  integer          :: i, j, ncoeffs

  write(6,'(/'' reading MO coefficients in AO basis'')')
  read(idmp,*)
  ncoeffs = 0

  do i = 1, nbsf
     do j = 1, nbsf

        read(idmp,'(es24.16)',err=10) vmoao(j,i)
        ncoeffs = ncoeffs + 1

     enddo
  enddo

  write(6,'(''  ====> complete'')')
  write(6,'(''  number of coefficients restored from file : '', &
         & i12)') ncoeffs

  return

10 continue
  write(6,'(/,'' ERROR::: something went wrong reading MO coeffs'')')
  write(6,'(''MO coeffs expected : '',i8,/,''MO coeff attempted to read : '',i8)') nbsf * nbsf, ncoeffs
  stop
  
end subroutine rdsamo

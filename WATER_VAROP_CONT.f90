subroutine GET_WATER_CONT_COEF(wr)
 use TELLURICA
 use OBSERVED
 use CONSTANTS, only: cl2
 use WATER_VAROP_CONT
 implicit none
 integer, intent(in) :: wr
 integer, parameter :: id_h2o = 1
 logical, parameter :: use_foreign_closure = .false.
 integer :: i, j, k
 integer, allocatable :: ileft(:)
 real(8), allocatable :: wnobs(:), frac(:)
 real(8), allocatable :: self_temp(:), foreign_ref(:)
 real(8) :: wn, dw, Tr, logTr
 real(8) :: xh2o, U, density_ratio
 real(8) :: selfc, foreignc, radiation_term

 if(allocated(water_cont_layer)) deallocate(water_cont_layer)
 if(allocated(water_cont_tau)) deallocate(water_cont_tau)

 allocate(water_cont_layer(nobs,nlayer), water_cont_tau(nobs))
 water_cont_layer = 0.d0; water_cont_tau = 0.d0
 
 allocate(wnobs(nobs), frac(nobs), ileft(nobs))
 allocate(self_temp(nwvcont), foreign_ref(nwvcont))
 dw = wavenumbers(2) - wavenumbers(1)

 if(use_foreign_closure) then
    foreign_ref = for_closure_absco_ref
 else
    foreign_ref = for_absco_ref
 endif
 do i = 1, nobs
    wn = 1.d8 / wobs(i)
    wnobs(i) = wn
    ! Do not silently extrapolate beyond the supplied table.
    if(wn < 0.d0 .or. wn > wavenumbers(nwvcont)) then
       write(6,'(a,i10,a,es20.10)') 'GET_WATER_CONT_COEF: point ', i, ' outside continuum table, wn = ', wn
       stop 'GET_WATER_CONT_COEF: unsupported spectral range'
    endif
    k = int((wn - wavenumbers(1)) / dw) + 1
    k = max(1, min(k, nwvcont - 1))
    ileft(i) = k
    frac(i) = (wn - wavenumbers(k)) / (wavenumbers(k+1) - wavenumbers(k))
 enddo
 do j = 1, nlayer
    xh2o = abudance(id_h2o,j)
    U = abudance_mol(id_h2o,j)
    if(U == 0.d0 .or. press(j) == 0.d0) cycle
    Tr = water_cont_ref_temp / teff(j)
    logTr = log(Tr)
    density_ratio = (press(j) / water_cont_ref_press_pa) * Tr
    self_temp = self_absco_ref * exp(self_texp * logTr)
    do i = 1, nobs
       k = ileft(i)
       wn = wnobs(i)
       selfc = (1.d0 - frac(i))*self_temp(k) + frac(i)*self_temp(k+1)
       foreignc = (1.d0 - frac(i))*foreign_ref(k) + frac(i)*foreign_ref(k+1)
       radiation_term = wn * tanh(0.5d0*cl2*wn/teff(j))
       water_cont_layer(i,j) = U * radiation_term * density_ratio * (xh2o*selfc + (1.d0-xh2o)*foreignc)
    enddo
 enddo
 water_cont_tau = sum(water_cont_layer, dim=2)
 deallocate(wnobs, frac, ileft)
 deallocate(self_temp, foreign_ref)

 if(wr == 1) write(6,'(a)') '  Preparation of water vapor continuum by layers completed successfully'

end subroutine GET_WATER_CONT_COEF
subroutine effective_depths(wavelength) ! Calculate effective depths of formation the residual in given spectrum point 
use radiation_field
use synth_data
use atmosphere_model, only : Nrhox
real(8) tau_eff_norm,d   ! Normalization factor for calculate average depth as tau5000 and depth's number
real(8) wavelength,kc,c10,kt  ! Wavelength (A),continuum opacity and log(10)
!real(8), allocatable :: kR(:),tauR(:),SR(:),y(:) ! "opacity", "optical depth","sources function" for residual intensity and array for numeric integration
real(8), allocatable :: depths_numbers(:)        ! atmosphere depths numbers for interpolation to found effective depth
integer(2), save :: first_time =0
real(8) expi
save 

  if(first_time == 0)then  ! Allocate arrays and calculate optical depths in continuum at 5000A 
   first_time=1
!   allocate (kR(Nrhox),SR(Nrhox),y(Nrhox),tauR(Nrhox),depths_numbers(Nrhox),stat=ios)
   allocate (depths_numbers(Nrhox),stat=ios)
   if(ios /= 0)stop 'Memory allocation for effective depths failed'
   c10=log(10.d0)
   do j=1,Nrhox
    depths_numbers(j)=j
   enddo
  endif

 call get_continuum_radiation_field(wavelength) ! Calculate contribution function in version Achmad et.al,1991, formulaes (10),(11)
 tau_eff_norm=0.d0; effective_depth_tau=0.d0; effective_depth_number=0.d0
!  do j=1,Nrhox                                 
!   kc=acont(j)+sigmac(j)                        ! total opacity in continuum
!   kr(j)=Aline(j)+4.d0*kc*S_cont(j)/H_cont(j)   ! R "opacity" from (10)
! R "source function"
!   Sr(j)=(1.d0/(kc*S_cont(j)+Aline(j)*H_cont(j)/4.d0))*(kc*((S_cont(j)-Snu(j))+(Jnu(j)-Hnu(j)/H_cont(j)*J_cont(j)))+Aline(j)*((H_cont(j)-Hnu(j)/4.d0+(Jnu(j)-Snu(j)))))
!   y(j)=c10*tau5000(j)*kr(j)/kappa5000(j)       ! integrand to R "optical depths"
!  enddo
!  call integ(logtau5000,y,tauR,Nrhox,logtau5000(1)*y(1)) ! R "optical depths"
 do j=1,Nrhox
!  CF=c10*tau5000(j)*kr(j)/kappa5000(j)*Sr(j)*exp(-tauR(j))  ! Contribution function from (20) notation
!  kc=acont(j)+sigmac(j)+AHline(j); kt=Aline(j)+kc 
  kc=acont(j)+sigmac(j); kt=Aline(j)+kc 
  CF=c10*tau5000(j)*(kt*Snu(j)-kc*S_cont(j))/kappa5000(j)*expi(2,taunu(j)) ! Contribution function from eq. (31)
  tau_eff_norm=tau_eff_norm+CF
!  effective_depth_tau=effective_depth_tau+CF*logtau5000(j)
  effective_depth_tau=effective_depth_tau+CF*tau5000(j)
!  write(403,*)j,cf; write(404,*)logtau5000(j),CF
 enddo

 if(tau_eff_norm <=0.d0.or.effective_depth_tau(1) <= 0.d0)then
   effective_depth_tau=0.d0
  else
   effective_depth_tau=log10(effective_depth_tau/tau_eff_norm)
  endif
   call linter(logtau5000,depths_numbers,Nrhox,effective_depth_tau,effective_depth_number,1)
 return
entry effective_depths_deallocation
! deallocate (kR,y,SR,tauR,depths_numbers)
 deallocate (depths_numbers)
end
 
 
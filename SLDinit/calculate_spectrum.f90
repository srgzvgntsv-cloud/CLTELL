subroutine calculate_spectrum ! Calculate synthetic spectrum using preliminary calculated opacities
use synth_data
use atmosphere_model
use radiation_field
use lines_data
use dflib
real(8) wavelength,Hnu_c,Inu_c(7),residual ! Wavelength(A),continuum flux and specific intensities; Hnu/Hnu_c
real(8) Inu_l(7)                           ! Array for save intensities in 7 rings (especially for Vrad /= 0)
real(8) Fc_first,Fc_mid,Fc_last,hnukt      ! continuum fluxes for prf output
real(8) shift,wave,flux,flux_c                         ! Vrad wavelengths shift and shifted wavelength, flux as Inu moment
real(8), allocatable :: specific_intensities(:,:),continuum_intensities(:,:),e_depths(:),n_depths(:)

 if(Vrad_switch > 0)then  ! Buffer for saving the specific intensities
  allocate (specific_intensities(7,number_of_wavelengths),continuum_intensities(7,number_of_wavelengths),&
            e_depths(number_of_wavelengths),n_depths(number_of_wavelengths),stat=i1)
  if(i1 /= 0)stop 'Not enough memory, decrease wavelength interval' 
 endif   

  allocate(synth_wl(number_of_wavelengths), synth_R(number_of_wavelengths))

 shift=0.d0; nmu=1; if(Vrad_switch > 0)nmu=7; ndepths=1; if(Vrad_switch == 1)ndepths=Nrhox

 do i=1,number_of_wavelengths  ! Cycle over wavelengths
  wavelength=wavelength_first+(i-1)*wavelength_step ! Given wavelength
  do mu=1,nmu
   Aline=0.d0
   do n=1,ndepths
    if(Vrad_switch > 0)shift=Vrad(n)/2.997925d5*mu_fix(mu)*wavelength ! Doppler shift at given mu and depth
	wave=wavelength-shift
    freq=2.997925d18/(wave)               ! Frequency
    freqlg=log(freq)                                  ! log frequency
     do k=1,Nrhox
      EhvkT(K)=exp(-freq*hkT(k))
      stim(k)=1.d0-EhvkT(k)
      Bnu(k)=1.47439d-27*freq**3*EhvkT(k)/stim(k)*1.d-20 ! Plank's function
     enddo
    if(output_format == 'out')then   
     hnukt=1.43868D8/wave/T(1) 
     units_factor=exp(50.7649141D0-5.D0*LOG(wave))/(exp(hnukt)-1.D0)/Bnu(1) ! scale intensities in mout output to [erg/(s*cm^2*rad^2*A)]
    endif
    call  get_continuum_flux(wave)                ! Interpolate early calculated continuum flux and intensities to the given wavelength
    Hnu_c=Hnu(1); Inu_c=Inu
    if(Vrad_switch > 0)continuum_intensities(1:7,i)=Inu_c(1:7) ! Save continuum intensities
    call get_continuum_opacity(wave)              ! Interpolate early calculated continuum opacities to the given wavelength
if(wave > 3600.and.wave < 4200)then
call kapp
call Josh
Hnu_c=Hnu(1)
endif
    call Hlinop(wave)                             ! Hydrogen lines opacity
    if(atomic_lines + molecular_lines > 0)then
	 if(Vrad_switch > 0)then
       m=(wave-wavelength_first)/wavelength_step+1; if(m < 1)cycle; if(m > number_of_wavelengths)cycle
     else
	   m=i
     endif
     read(200,rec=m)aline_buffer
    else 
	 aline_buffer=0.d0
    endif
    Aline(n)=aline_buffer(n)+AHline(n)
   enddo
   if(ndepths == 1)Aline=aline_buffer+AHline
    call Josh; Inu_l(mu)=Inu(mu)
    if(mu == 1)then
	 call effective_depths(wave) ! calculate an effective stellar atmosphere depths given spectrum point are formed
     if(Vrad_switch /= 0)then
	  e_depths(i)=effective_depth_tau(1); n_depths(i)=effective_depth_number(1)
     endif
    endif
  enddo
  if(Vrad_switch > 0)specific_intensities(1:7,i)=Inu_l(1:7) ! Save specific intensities
  if(Vrad_switch == 0)then
  residual=Hnu(1)/Hnu_c
   select case (output_format)  ! Output calculated spectrum in different formats
   case ('my ')
    write(16,'(f10.3,f8.3,1p2e11.4,0pf12.3,f10.2)')wavelength,residual,Hnu(1),Hnu_c,effective_depth_tau,effective_depth_number
	write(26,rec=i+1)Hnu_c,Hnu(1),Inu_c,Inu,effective_depth_tau,effective_depth_number
   case ('prf')
    write(16,'(2f16.8)')wavelength,residual
   case ('out')
    write(16,'(8f12.5)')wavelength,Inu*units_factor
   end select
  endif

	synth_wl(i) = wavelength
	synth_R(i) = residual
 enddo


 if(Vrad_switch /= 0)then
  call disk_integration(specific_intensities,continuum_intensities,e_depths,n_depths,number_of_wavelengths)
  deallocate (specific_intensities,continuum_intensities,e_depths,n_depths)
 endif

 call effective_depths_deallocation

 if(output_format == 'prf')then     ! Especially for prf ouput - print continuum fluxes
  write(16,"(' Temperature=',f11.2,'     ,  Gravity=',f11.6)")Teff,glog
  call get_continuum_flux(wavelength_first); Fc_first=Hnu(1)
  call get_continuum_flux(wavelength_last); Fc_last=Hnu(1)
  wavelength=(wavelength_first+wavelength_last)/2.d0
  call get_continuum_flux(wavelength); Fc_mid=Hnu(1)
  write(16,"(' Fcblue=',e12.6,' Fcred=',e12.6,' Fcmid=',e12.6)")Fc_first,Fc_last,Fc_mid
 endif

end

subroutine disk_integration(specific_intensities,continuum_intensities,e_depths,n_depths,number_of_wavelengths) ! Calculate second moment from specific intensity."Smooth" spctrum
use radiation_field, only : mu_fix,effective_depth_tau,effective_depth_number,Bnu
use synth_data, only : Vrad,output_format,wavelength_first,wavelength_step
use atmosphere_model, only : T
implicit real(8) (a-h,o-z)
real(8) specific_intensities(7,number_of_wavelengths),continuum_intensities(7,number_of_wavelengths)
integer, parameter :: number_knots=50
real(8) intensity(number_knots),interpolated(number_knots),x(number_knots),weights(number_knots),xold(7),yold(7),shift,shift2,flux,flux_c
real(8) e_depths(number_of_wavelengths),n_depths(number_of_wavelengths)

 call MakeGaussMatrix(0.d0,1.d0,x,weights,number_knots) ! Calculate Gauss quadrature data x: knots, weights - quadrature weights

 do i=1,number_of_wavelengths
  wavelength=wavelength_first+(i-1)*wavelength_step
   do mu=1,number_knots
    shift=Vrad(1)/2.997925d5*x(mu)*wavelength
	do l=1,7
     shift2=Vrad(1)/2.997925d5*mu_fix(7-(l-1))*wavelength
     m=(wavelength-(shift-shift2)-wavelength_first)/wavelength_step+1; if(m < 1)m=1; if(m > number_of_wavelengths)cycle
     xold(l)=mu_fix(7-(l-1))
     yold(l)=specific_intensities(7-(l-1),m)
	enddo 
    call linter(xold,yold,7,x,interpolated,number_knots) ! Interpolate specific intensities to new mu grid
    intensity(mu)=interpolated(mu)
   enddo
  flux=0.d0
   do mu=1,number_knots
    flux=flux+x(mu)*intensity(mu)*weights(mu)
   enddo
   do mu=1,number_knots
    shift=Vrad(1)/2.997925d5*x(mu)*wavelength
     do l=1,7
      shift2=Vrad(1)/2.997925d5*mu_fix(7-(l-1))*wavelength
      m=(wavelength-(shift-shift2)-wavelength_first)/wavelength_step+1; if(m < 1)m=1; if(m > number_of_wavelengths)cycle
      xold(l)=mu_fix(7-(l-1))
      yold(l)=continuum_intensities(7-(l-1),m)
     enddo
    call linter(xold,yold,7,x,interpolated,number_knots)
    intensity(mu)=interpolated(mu)
   enddo
  flux_c=0.d0
   do mu=1,number_knots
    flux_c=flux_c+x(mu)*intensity(mu)*weights(mu)
   enddo
  flux=flux/4.d0; flux_c=flux_c/4.d0
  residual=flux/flux_c
   select case (output_format)  ! Output calculated spectrum in different formats
   case ('my ')
    write(16,'(f10.3,f8.3,1p2e11.4,0pf12.3,f10.2)')wavelength,residual,flux,flux_c,e_depths(i),n_depths(i)
	write(26,rec=i+1)flux_c,flux,continuum_intensities(1:7,i),specific_intensities(1:7,i),e_depths(i),n_depths(i)
   case ('prf')
    write(16,'(2f16.8)')wavelength,residual
   case ('out')
    hnukt=1.43868D8/wavelength/T(1) 
    units_factor=exp(50.7649141D0-5.D0*LOG(wavelength))/(exp(hnukt)-1.D0)/Bnu(1)
    write(16,'(8f12.5)')wavelength,specific_intensities(1:7,i)*units_factor
   end select
 enddo

end
 

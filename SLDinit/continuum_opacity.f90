subroutine continuum_opacity ! calculate continuum opacity in fixed wavelengths and interpolate it's to specteal wavelengths points
use synth_data
use continuum_wavelengths
use atmosphere_model
use radiation_field
real(8),allocatable ::  x_cont(:),y_cont(:) ! Arrays for interpolation: x_cont - wavelngths, y_cont - opacity
real(8),allocatable ::  continuum_opacities(:,:,:) ! First index: 1- absorbtion, 2- scattering,3-source function
real(8),allocatable ::  continuum_flux(:),continuum_intensities(:,:)
real(8), allocatable :: Snu_c(:,:),Hnu_c(:,:),Jnu_c(:,:) ! Save continuum radiation field data for calculation contribution function
real(8) wavelength
save



 number_continuum_points=2     ! Determine number of fixed continuum points inside calculation wavelengths interval
 do i=1,number_fixed_continuum_points
  if(wavelengths_continuum(i) < wavelength_first)cycle
  if(abs(wavelengths_continuum(i)-wavelength_first) < 1.d-3)cycle ! Strange, but <= doesn't work...
  if(wavelengths_continuum(i) > wavelength_last)exit
  if(abs(wavelengths_continuum(i)-wavelength_last) < 1.d-3)exit
  number_continuum_points=number_continuum_points+1
 enddo

 if(number_continuum_points == 2)number_continuum_points=3
   allocate(x_cont(number_continuum_points),y_cont(number_continuum_points),continuum_opacities(3,number_continuum_points,Nrhox),&
           continuum_flux(number_continuum_points),continuum_intensities(7,number_continuum_points),&
		   Snu_c(Nrhox,number_continuum_points),Hnu_c(Nrhox,number_continuum_points),Jnu_c(Nrhox,number_continuum_points),stat=ios)
  if(ios /= 0)stop 'Memory allocation in continuum opacity failed'
  x_cont(1)=wavelength_first; x_cont(number_continuum_points)=wavelength_last
   if(number_continuum_points == 3)then
    x_cont(2)=(wavelength_first+wavelength_last)/2.d0
   else
    m=1
    do i=1,number_fixed_continuum_points-1
     if(wavelengths_continuum(i) < wavelength_first)cycle
     if(abs(wavelengths_continuum(i)-wavelength_first) < 1.d-3)cycle 
	 if(wavelengths_continuum(i) > wavelength_last)exit
     if(abs(wavelengths_continuum(i)-wavelength_last) < 1.d-3)exit
	 m=m+1
     x_cont(m)=wavelengths_continuum(i)
    enddo
   endif
 
 do i=1,number_continuum_points ! Calculate continuum opacity in selected wavelengths points
  freq=2.997925d18/x_cont(i)
  freqlg=log(freq)
   do k=1,Nrhox
    EhvkT(K)=exp(-freq*hkT(k))
    stim(k)=1.d0-EhvkT(K)
    Bnu(k)=1.47439d-27*freq**3*EhvkT(k)/stim(k)*1.d-20
   enddo 
   call Kapp
   continuum_opacities(1,i,1:Nrhox)=Acont(1:Nrhox)
   continuum_opacities(2,i,1:Nrhox)=sigmaC(1:Nrhox)
   continuum_opacities(3,i,1:Nrhox)=Scont(1:Nrhox) ! This one if for NLTE synthetic spectra
! Calculate continuum output flux in fixed wavelengths points
   Aline=0.d0 ! Set lines absorbtion equial zero
   call Josh
   continuum_flux(i)=Hnu(1)
   continuum_intensities(1:7,i)=Inu(1:7)
   Snu_c(1:Nrhox,i)=Snu(1:Nrhox);Hnu_c(1:Nrhox,i)=Hnu(1:Nrhox);Jnu_c(1:Nrhox,i)=Jnu(1:Nrhox);
 enddo
return

entry get_continuum_opacity(wavelength) ! Interpolate continuum opacity to given wavelength
 
 do j=1,Nrhox
  y_cont(1:number_continuum_points)=continuum_opacities(1,1:number_continuum_points,j)
  id=map1(x_cont,y_cont,number_continuum_points,wavelength,Acont(j),1)
  y_cont(1:number_continuum_points)=continuum_opacities(2,1:number_continuum_points,j)
    y_cont(1:number_continuum_points)=continuum_opacities(2,1:number_continuum_points,j)
  id=map1(x_cont,y_cont,number_continuum_points,wavelength,sigmaC(j),1)
  y_cont(1:number_continuum_points)=continuum_opacities(3,1:number_continuum_points,j)
  id=map1(x_cont,y_cont,number_continuum_points,wavelength,Scont(j),1)
 enddo
return

entry get_continuum_radiation_field(wavelength)
 
 do j=1,Nrhox 
  y_cont(1:number_continuum_points)=Snu_c(j,1:number_continuum_points)
   id=map1(x_cont,y_cont,number_continuum_points,wavelength,S_cont(j),1)
  y_cont(1:number_continuum_points)=Hnu_c(j,1:number_continuum_points)
   id=map1(x_cont,y_cont,number_continuum_points,wavelength,H_cont(j),1)
  y_cont(1:number_continuum_points)=Jnu_c(j,1:number_continuum_points)
   id=map1(x_cont,y_cont,number_continuum_points,wavelength,J_cont(j),1)
 enddo
return

entry get_continuum_flux(wavelength) ! Interpolate continuum opacity to given wavelength
 
  id=map1(x_cont,continuum_flux,number_continuum_points,wavelength,Hnu(1),1)
  do i=1,7
   y_cont(1:number_continuum_points)=continuum_intensities(i,1:number_continuum_points) 
   id=map1(x_cont,y_cont,number_continuum_points,wavelength,Inu(i),1)
 enddo
return

entry print_continuum_intensities ! especially entry to mout output
  write(16,'(i12)')number_continuum_points
   do i=1,number_continuum_points
    freq=2.997925d18/x_cont(i)
    EhvkT(1)=exp(-freq*hkT(1))
    stim(1)=1.d0-EhvkT(1)
    Bnu(1)=1.47439d-27*freq**3*EhvkT(1)/stim(1)*1.d-20 
    hnukt=1.43868D8/x_cont(i)/T(1) 
    units_factor=exp(50.7649141D0-5.D0*LOG(x_cont(i)))/(exp(hnukt)-1.D0)/Bnu(1) ! scale intensities in mout output to [erg/(s*cm^2*rad^2*A)]
    write(16,'(8f12.5)')x_cont(i),(continuum_intensities(k,i)*units_factor,k=1,7)
!    write(16,'(f12.5,1p7e12.5)')x_cont(i),(continuum_intensities(k,i),k=1,7)
   enddo
return
   	 
entry deallocate_continuum_arrays
 deallocate(x_cont,y_cont,continuum_opacities,continuum_flux,continuum_intensities,Snu_c,Hnu_c,Jnu_c)

end

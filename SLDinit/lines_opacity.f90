subroutine lines_opacity ! Preselect spectral lines and calculate opacity in atomic and (optionaly) molecular lines
use synth_data
use lines_data
use chemical_elements
use dflib
use radiation_field, only : mu_fix,Inu
use mol_dat, only : Mol_id, Nmol
integer integer_wavelength,integer_E ! Data in binary file are placed as a integer(4) values
integer(2) integer_gf,integer_code,integer_dumping1,integer_dumping2,integer_dumping3 ! Data in binary file are placed as a integer(2) values
real(8) E,limb_darkening,a,b,wavelength
character(5) new_VALD
character(100)string

  if(atomic_lines == 1)then ! Take atomic lines data from binary file
	nrec=first_record_atomic-1; mode=1
	do
		read(3,rec=nrec)integer_wavelength,integer_gf,integer_E,integer_code,integer_dumping1,integer_dumping2,integer_dumping3
		central_wavelength=integer_wavelength/1000.d0 ! (A)
		if(central_wavelength > wavelength_last+2.5d0)exit ! Complete lines selection
		nrec=nrec+1
		gf=10**(integer_gf/1000.d0); if(gf > 5)gf=-gf ! For prevent possible bugs in source VALD list
		code=integer_code/100.d0; E=integer_E/1000.d0; Elow=E*2.997925d10 !Convert integers to float (true) values
		gamma_Rad=10**(integer_dumping1/100.d0) ;gamma_Stark=10**(integer_dumping2/100.d0) ;gamma_WdW=10**(integer_dumping3/100.d0)
		atomic_number=code+0.5; ionization=integer_code-atomic_number*100+1
		write(line_identifier,'(a2,i2)')Elem_id(atomic_number),ionization
		if(ionization > 6)cycle ! Maximum 5th ions
		call line_absorbtion(mode) ! Calculate opacity, due by line in all wavelengths. Mode=1 mean that line is atomic
	enddo
! Type first record in prf or out file
	if(output_format /= 'my ')then
		rewind (16)
		write(300, "(i4,'- number of spectral lines')")number_of_lines  
		do i=1,number_of_lines
			read(16,'((a))')string
			write(300,'((a))')trim(string)
		enddo
		rewind(300); rewind(16)
		do i=1,number_of_lines+1
			read(300,'((a))')string
			write(16,'((a))')trim(string)
		enddo
		close (300,disp='delete')
	endif
  endif

  if(atomic_lines == 2)then ! Take atomic lines data from VALD file
   read(3,*,iostat=ios)wstart,wmax,Nlines
   if(ios /= 0)stop 'Wrong data in VALD file'
   read(3,*) ! Skip head's records
   read(3,*)
   do nl=1,Nlines ! Read lines data
    mode=1; new_VALD=' '; VALD_string(1:200)=' '
	read(3,'((a))',iostat=ios)VALD_string
	if(ios /= 0)stop 'Number of lines in VALD file is incorrect'
	read(VALD_string,*,iostat=ios)new_VALD,central_wavelength,E,dummy,gf,gamma_Rad,gamma_Stark,gamma_WdW
	if(ios /= 0)then
	 write(*,*)'Wrong record number ',nl,' in VALD file'; stop
    endif
	if(new_VALD == 'H  1')cycle  ! Skip Hydrogen lines
	if(central_wavelength < wavelength_first-2.5d0)cycle ! skip lines outside wavelengths region
    if(central_wavelength > wavelength_last+2.5d0)exit
	if(if_ev == 0)E=E/8066.d0
    Elow=E*2.997925d10*8066.d0
    gf=10**gf
	line_identifier=new_VALD(1:4)
	if(molecular_lines == 2)then
     do l=1,Nmol
	  k=index(new_VALD,trim(Mol_id(l)))
	  if(k /= 0)then
       mode = 2; exit
      endif
     enddo
	endif
	if(mode == 2)then
	 if(new_VALD(3:3) == ' '.and.new_VALD(4:4) == '1')line_identifier(4:4)=' '  ! CO, C2 etc.
	 if(new_VALD(3:3) == ' '.and.new_VALD(4:4) == '2')line_identifier(3:4)='+ ' ! CO+,CH+ etc.
	 if(new_VALD(3:3) /= ' '.and.new_VALD(5:5) == '2')line_identifier(4:4)='+'  ! MgH+ etc.
	 if(gamma_Rad == 0.d0)then
	  gamma_Rad=2.223d13/(central_wavelength**2)
     else
	  gamma_Rad=10**gamma_Rad
     endif
     gamma_Stark=1.d-5   ! Stark and Wan der Waals constant fixed for molecular lines
     gamma_WdW=1.d-7
	else  
     read(line_identifier(3:4),*)ionization
	 if(ionization > 6)cycle ! Maximum 5th ions
     do i=1,99  ! Search the atomic number
      if(line_identifier(1:2) == Elem_id(i))exit
     enddo
     atomic_number=i
	 if(atomic_number == 1)cycle  ! Skip Hydrogen lines  
     code=atomic_number+(ionization-1)/100.d0
  	 call calculate_line_constants(central_wavelength,E,code,gamma_Rad,gamma_Stark,gamma_WdW) ! Calculate dumping constants if in VALD file it's equial zero
     gamma_Rad=10**gamma_Rad; gamma_Stark=10**gamma_Stark; if(gamma_WdW < 10)gamma_WdW=10**gamma_WdW
	endif
    call line_absorbtion(mode) ! Calculate opacity, due by line in all wavelengths
   enddo
   if(molecular_lines == 2.and.output_format /= 'my ')then  ! Type first record in prf or out file
    rewind (16)
    write(300, "(i4,'- number of spectral lines')")number_of_lines  
    do i=1,number_of_lines
	 read(16,'((a))')string
	 write(300,'((a))')trim(string)
    enddo
	rewind(300); rewind(16)
	do i=1,number_of_lines+1
	 read(300,'((a))')string
	 write(16,'((a))')trim(string)
    enddo
	close (300,disp='delete')
   endif
 
  endif

  if(molecular_lines == 1)then  ! Molecular lines from binary format file
   nrec=first_record_molecular-1; mode=2
    do 
     read(13,rec=nrec,iostat=ios)line_identifier,integer_wavelength,integer_gf,integer_E,integer_dumping1
     if(ios /= 0)exit
     central_wavelength=integer_wavelength/10000.d0  ! (A)
     if(central_wavelength > wavelength_last+2.5d0)exit
     nrec=nrec+1
     gf=10**(integer_gf/1000.d0); E=integer_E/1000.d0; Elow=E*2.997925d10 !Convert integers to float (true) values
     gamma_Stark=1.d-5   ! Stark and Wan der Waals constant fixed for molecular lines
     gamma_WdW=1.d-7
     call line_absorbtion(mode) ! Calculate opacity, due by line in all wavelengths. Mode=2 mean that line is molecular
    enddo
	call sort_id_lines
   endif

    if(output_format == 'out')then ! out output
     write(16,*)'           7'
	 write(16,'(f24.5,6f12.5)')mu_fix
	 call print_continuum_intensities
	 write(16,'(i12)')number_of_wavelengths
	else                            ! prf output
	 if(output_format == 'prf')then
	  wavelength=(wavelength_first + wavelength_last)/2.d0 ! Calculate continuum limb darkening for middle wavelength
	  call get_continuum_flux(wavelength)
      a=1.d0-mu_fix(4)
      b=1.d0-mu_fix(7)
      limb_darkening=(a*(1.d0-Inu(4)/Inu(1))+b*(1.d0-Inu(7)/Inu(1)))/(a*a+b*b) ! Interpolation via disk rings 
      write(16,"(i6,f8.4,'- number of wavelegths,limb darkening')")number_of_wavelengths,limb_darkening
     endif
    endif

end

subroutine line_absorbtion(mode) ! Calculate absorbtion due given spectral line. mode different for atomic and molecular lines
use synth_data
use atmosphere_model
use radiation_field
use lines_data
use temporary_arrays
use Hetab
use chemical_elements
use mol_dat
use HlineTab, only : wAir,nHydSeries,nHydLines
implicit real(8) (a-h,o-z)
real(8) id_Hydrogen(80)  ! wavelengths of n_id_Hydrogen lines inside wavelengths of computed spectrum, which will be putted in output lines list
real(8) number_density,mass
integer(2) id_Hydrogen_l(80)
integer(2), save :: first_time =0
integer(2) Stark_switch
character(5) new_VALD
save

 if(first_time == 0)then  ! Precalculate values, need for lines opacity calclulations
  first_time=1
   do j=1,Nrhox
    Tx(j)=(xDens2(j,1,1)+0.42d0*rho(j)*1.77245d0/0.026538d0*&
           xDens2(j,2,1))*(T(j)/10000.d0)**0.3d0
    vt2(j)=Vturb(j)**2
   enddo 

   do i=1,number_of_wavelengths ! Initialize temporary file, contained total lines absorbtion
    Aline=0.d0
	 write(200,rec=i)Aline
   enddo

   n_id_Hydrogen=0; id_Hydrogen_l=0 ! save Hydrogen lines for put it's in output lines list
   do i=1,nHydSeries
    do j=1,nHydLines
	 if(wAir(j,i) >= wavelength_first.and.wAir(j,i) <= wavelength_last)then
	  n_id_Hydrogen=n_id_Hydrogen+1; id_Hydrogen(n_id_Hydrogen)=wAir(j,i); id_Hydrogen_l(n_id_Hydrogen)=1
     endif
    enddo
   enddo
 endif
!-----------------------------------------------------------------------
  Stark_switch=0
  if(atomic_number == 2.and.ionization == 1.and.Teff > 9000)then   !  Especialy for 20 HeI lines. In future extend list to all lines with specific Stark dumping
   iw=central_wavelength
   do i=1,23
    if(iw == wl_He(i))then
	  Stark_switch=i;  if(i == 21)Stark_switch=20; exit
    endif
   enddo
  endif

  fr=2.997925d18/central_wavelength ! Frequency
  fc=fr/2.997925d10
  wn=1.d-7*fc*fr

  if(mode == 1)then
   mass=Atmass(atomic_number)
  else
   nm=0
   do n=1,Nmol
    if(line_identifier == Mol_id(n))then
	 nm=n; exit
	endif
   enddo 
   if(nm == 0)return ! Such molecule not included
   mass=Mol_mass(nm)
  endif

  if(Stark_switch /= 0)gamma_Stark=0.d0

   do j=1,Nrhox

    TDop=sqrt(Tk(j)*1.2048192d14/mass+vt2(j))*1.d5

    Dopnu=1.d0/(TDop*fc)
    if(mode == 1)then
	 number_density=xDens(j,atomic_number,ionization)
    else
	 number_density=x_mol(j,nm)
    endif
    al_int(j)=exp(-Elow*HkT(j))*(1.d0-exp(-fr*HkT(j)))*number_density*Dopnu*gf ! Integral line absorbtion
	if(gamma_WdW < 10)then
	 gWdw=gamma_WdW*Tx(j)
    else  ! Anstee & O'Mara 
	 sigma=int(gamma_WdW)*2.8002856d-21 ! (a0=5.29177249d-11)**2
	 pWdW=gamma_WdW-int(gamma_WdW)
     G_X=(4.d0-pWdW)/2.d0-1.d0
     GAMMAF=1.d0+(-0.5748646d0+(0.9512363d0+(-0.6998588d0+(0.4245549d0-0.1010678d0*G_X)*G_X)*G_X)*G_X)*G_X
	 v0=1.d4  ! 10km/sec
	 gWdw=(4.d0/3.1415925)**(pWdW/2.d0)*GAMMAF*v0*sigma    ! Gray, p.247. Check - lok like VBAR instead v0!
     vbar=sqrt(T(j)*2.117273132059405d4*(1.d0/1.008d0+1.d0/mass))
	 gWdW=gWdW*(VBAR/v0)**(1.d0-pWdW)
	 gWdW=gWdW*(xDens2(j,1,1)+0.42d0*rho(j)*1.77245d0/0.026538d0*xDens2(j,2,1))*2.d0*1.d6
	endif
    Adump(j)=(gamma_Rad+gamma_Stark*xNe(j)+gWdW)*0.0795773d0*Dopnu  ! Pressure dumping

    Dopdump(j)=Dopnu*wn                                                        ! Doppler dumping 
   enddo

  call get_continuum_opacity(central_wavelength)
  Abtot=Acont+Sigmac ! Total continuum opacity
  w=Abtot

  call line_absorbtion_value(0,Stark_switch,if_stop) ! Intensity in centre of given line
  if(if_stop == 1)return
   call  get_continuum_flux(central_wavelength)                ! Interpolate early calculated continuum flux and intensities to the given wavelength
   r=Hnu(1)
  freq=2.997925d18/central_wavelength
  freqlg=log(freq)
   do k=1,Nrhox
    EhvkT(K)=exp(-freq*hkT(k))
    stim(k)=1.d0-EhvkT(K)
    Bnu(k)=1.47439d-27*freq**3*EhvkT(k)/stim(k)*1.d-20
   enddo 

 !  call Hlinop(central_wavelength)                             ! Hydrogen lines opacity
   Aline=aline_buffer !+AHline  !Not need for identification: invisible lines - H lines wings   ! given line opacity
   call Josh

   r=1.d0-Hnu(1)/r                                             ! line centre depth
   rcentre=1.d0-r
   if(r > 0.0002d0)then                                            ! Skip identification output for invisible lines
   if_Hydrogen=0
   if(n_id_Hydrogen /= 0)then
    do i=1,n_id_Hydrogen
	 if(central_wavelength > id_Hydrogen(i).and.id_Hydrogen_l(i) == 1)then
	  if_Hydrogen=i; id_Hydrogen_l(i)=0; exit
	 endif
    enddo
   endif  
                                              
    if(output_format == 'my ')then
	 if(if_Hydrogen /= 0)then
	  write(27,'(f10.4,1x,a4,f8.4)')id_Hydrogen(if_Hydrogen),'H  1',1.d0-r
      number_of_lines=number_of_lines+1
     endif
     write(27,'(f10.4,1x,a4,f8.4)')central_wavelength,line_identifier,1.d0-r ! Lines identification in .lin file
    else                                                                     ! same in stupid prf or out format file
     g_Lande=99.000d0; if(atomic_lines == 2)read(VALD_string,*)new_VALD,d,d,d,d,d,d,d,g_Lande
	 if(g_lande == 0.d0)g_Lande=99.000d0
	 if(if_Hydrogen /= 0)then
	  write(16,"(f10.4,2h,',a4,2h',,f8.4,',',f4.1,',',f7.3,',', F8.3,',',F8.3,',',F8.3,',',F8.3,',',F7.3,', ',3h' ')")&
      id_Hydrogen(if_Hydrogen),'H  1',Elow/2.997925d10/8066.d0,V_micro,log10(gf),log10(gamma_Rad),log10(gamma_Stark),log10(gamma_WdW),g_Lande,r
      number_of_lines=number_of_lines+1
     endif
	 if(gamma_Stark == 0.d0)gamma_Stark=1.d0  ! Prevent crash for non-standard Stark dumpening (log10(0))
	 if(gamma_WdW < 10)gamma_WdW=log10(gamma_WdW)
     write(16,"(f10.4,2h,',a4,2h',,f8.4,',',f4.1,',',f7.3,',', F8.3,',',F8.3,',',F8.3,',',F8.3,',',F7.3,', ',3h' ')")&
     central_wavelength,line_identifier,Elow/2.997925d10/8066.d0,V_micro,log10(gf),log10(gamma_Rad),log10(gamma_Stark),gamma_WdW,g_Lande,r
    endif
    number_of_lines=number_of_lines+1
   endif
  
  if_stop=0
  index_wavelength=(central_wavelength-wavelength_first)/wavelength_step+1 ! record number in temporary absorbtion file
  if(index_wavelength < 1)index_wavelength = 1
  if(index_wavelength > number_of_wavelengths)index_wavelength = number_of_wavelengths
  call line_absorbtion_value(index_wavelength,Stark_switch,if_stop) ! Absorbtion in centre of given line
  if(if_stop == 1)return

  if(index_wavelength > 1)then  ! If happened given lines due only in wavelengths limit
   do i=index_wavelength-1,1,-1  ! Calculate absorbtion profile at blue wing
    call line_absorbtion_value(i,Stark_switch,if_stop)
    if(if_stop /= 0)then ! end of given line profile
	 if(if_database /= 0)ww1= wavelength_first+(i-1)*wavelength_step
	 exit
    endif
   enddo
  endif

   if(index_wavelength < number_of_wavelengths-1)then  ! If happened given lines due only in wavelengths limit
   do i=index_wavelength+1,number_of_wavelengths  ! Calculate absorbtion profile at red wing
    call line_absorbtion_value(i,Stark_switch,if_stop)
    if(if_stop /= 0)then
	 if(if_database /= 0)then
	  ww2=wavelength_first+(i-1)*wavelength_step
	  if(rcentre < 0.9998d0.and.ww1 /= 0.d0.and.ww2 /= 0.d0)then
	   call line_profile(ww1,ww2,central_wavelength,line_identifier,rcentre,Stark_switch)
      endif
     endif
	 exit
	endif
   enddo
  endif
end

subroutine line_absorbtion_value(index_wavelength,Stark_switch,if_stop) ! Calculate line opacity in given wavelength
use synth_data
use lines_data
use radiation_field
use atmosphere_model
use Hetab
use temporary_arrays, only : w
real(8) wavelength,x,max_ratio,kappa,xDop     ! wavelength,distance from line centre,ratio to continuum opacity,opacity, Doppler shift
integer(2) Stark_switch
real(8) Voigt

  if_stop=0
  if(index_wavelength == 0)then
   wavelength=central_wavelength
  else
   wavelength=wavelength_first+(index_wavelength-1)*wavelength_step
  endif
  if(wavelength < wavelength_first-2.5d0.or.wavelength > wavelength_last+2.5d0)then
   if_stop=1;return
  endif

  max_ratio=-1.d30                         ! Maximum of line/continuum absorbtion ratio
   do j=1,Nrhox
    x=wavelength-central_wavelength             ! Distance from line centre
    if(Stark_switch == 0)then
     xDop=abs(x)/10.d0*Dopdump(j) 
	 kappa=al_int(j)*Voigt(xDop,Adump(j))
    else
	 call He_prof(Stark_switch,j,x,kappa)
    endif
	aline_buffer(j)=kappa                              ! Given line centre absorbtion vs. depth
     ratio=kappa/w(j)
	 if(ratio > max_ratio)max_ratio=ratio
    enddo
   if(max_ratio < 0.001d0)then  ! Skip lines, which absorbtion is always smaller then 1% from continuum opacity
    if_stop=1; return
   endif
   if(index_wavelength == 0)return
    
   read(200,rec=index_wavelength)Aline; Aline=Aline+aline_buffer
   write(200,rec=index_wavelength)Aline

end

subroutine sort_id_lines ! Sort atomic and molecular lines by wavelengths
use synth_data
character(88), allocatable :: string(:)
character(88) str
real(8), allocatable :: wavelength(:)
real(8) wbuf

 allocate (string(number_of_lines),wavelength(number_of_lines))
 lunit=27; if(output_format /= 'my ')lunit=16; 
 rewind (lunit)
 do i=1,number_of_lines    ! Find first record with molecular line
  read(lunit,'((a))')string(i)
  read(string(i),*)wavelength(i)
 enddo
 rewind (lunit)
 do i=number_of_lines-1,1,-1  ! Sort by wavelengths
  do j=1,i
   if(wavelength(j) > wavelength(j+1))then
    wbuf=wavelength(j)
    str=string(j)
    wavelength(j)=wavelength(j+1); wavelength(j+1)=wbuf
    string(j)=string(j+1); string(j+1)=str
   endif
  enddo
 enddo 
 if(output_format /= 'my ')write(16, "(i4,'- number of spectral lines')")number_of_lines ! write number of lines
  do i=1,number_of_lines
   write(lunit,'((a))')trim(string(i))
  enddo
 
  deallocate(string,wavelength)

end


 subroutine calculate_line_constants(w,e,code,gr,gs,gw)
 implicit real(8) (a-h,o-z)
  DIMENSION PotIon(594)
  DATA (PotIon(i),i=1,114)    /&
   13.598,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,&               ! 1
   24.587,  54.416,   0.   ,   0.   ,   0.   ,   0.   ,&               ! 2
   5.392,  75.638, 122.451,   0.   ,   0.   ,   0.    ,&               ! 3
   9.322,  18.211, 153.893, 217.713,   0.   ,   0.    ,&               ! 4
   8.298,  25.154,  37.930, 259.368, 340.217,   0.    ,&               ! 5
   11.260,  24.383,  47.887,  64.492, 392.077, 489.981,&               ! 6
   14.534,  29.601,  47.448,  77.472,  97.888, 552.057,&               ! 7
   13.618,  35.116,  54.934,  77.412, 113.896, 138.116,&               ! 8
   17.422,  34.970,  62.707,  87.138, 114.240, 157.161,&               ! 9
   21.564,  40.962,  63.45 ,  97.11 , 126.21 , 157.93 ,&               !10
    5.139,  47.286,  71.64 ,  98.91 , 138.39 , 172.15 ,&               !11
    7.646,  15.035,  80.143, 109.24 , 141.26 , 186.50 ,&               !12
    5.986,  18.828,  28.447, 119.99 , 153.71 , 190.47 ,&               !13
    8.151,  16.345,  33.492,  45.141, 166.77 , 205.05 ,&               !14
   10.486,  19.725,  30.18 ,  51.37 ,  65.023, 220.43 ,&               !15
   10.36 ,  23.33 ,  34.83 ,  47.30 ,  72.68 ,  88.049,&               !16
   12.967,  23.81 ,  39.61 ,  53.46 ,  67.80 ,  97.03 ,&               !17
   15.759,  27.629,  40.74 ,  59.81 ,  75.02 ,  91.007,&               !18
    4.341,  31.625,  45.72 ,  60.91 ,  82.66 , 100.0  /                !19
  DATA (PotIon(i),i=115,228)    /&
    6.113,  11.871,  50.908,  67.10 ,  84.41 , 108.78 ,&               !20
    6.54 ,  12.80 ,  24.76 ,  73.47 ,  91.66 , 111.1  ,&               !21
    6.82 ,  13.58 ,  27.491,  43.266,  99.22 , 119.36 ,&               !22
    6.74 ,  14.65 ,  29.310,  46.707,  65.23 , 128.12 ,&               !23
    6.766,  16.50 ,  30.96 ,  49.1  ,  69.3  ,  90.56 ,&               !24
    7.435,  15.640,  33.667,  51.2  ,  72.4  ,  95.   ,&               !25
    7.870,  16.18 ,  30.651,  54.8  ,  75.0  ,  99.0  ,&               !26
    7.86 ,  17.06 ,  33.50 ,  51.3  ,  79.5  , 102.0  ,&               !27
    7.635,  18.168,  35.17 ,  54.9  ,  75.5  , 108.0  ,&               !28
    7.726,  20.292,  36.83 ,  55.2  ,  79.9  , 103.0  ,&               !29
    9.394,  17.964,  39.722,  59.4  ,  82.6  , 108.0  ,&               !30
    5.999,  20.51 ,  30.71 ,  64.0  ,   0.   ,   0.   ,&               !31
    7.899,  15.934,  34.22 ,  45.71 ,  93.5  ,   0.   ,&               !32
    9.81 ,  18.633,  28.351,  50.13 ,  62.63 , 127.6  ,&               !33
    9.752,  21.19 ,  30.820,  42.944,  68.3  ,  81.70 ,&               !34
   11.814,  21.8  ,  36.   ,  47.3  ,  59.7  ,  88.6  ,&               !35
   13.999,  24.359,  36.95 ,  52.5  ,  64.7  ,  78.5  ,&               !36
    4.177,  27.28 ,  40.0  ,  52.6  ,  71.0  ,  84.4  ,&               !37
    5.695,  11.030,  43.6  ,  57.0  ,  71.6  ,  90.8  /                !38
  DATA (PotIon(i),i=229,342)    /&
    6.38 ,  12.24 ,  20.52 ,  61.8  ,  77.0  ,  93.0  ,&               !39
    6.84 ,  13.13 ,  22.99 ,  34.34 ,  81.50 ,   0.   ,&               !40
    6.88 ,  14.32 ,  25.04 ,  38.3  ,  50.55 , 102.6  ,&               !41
    7.099,  16.15 ,  27.16 ,  46.4  ,  61.2  ,  68.0  ,&               !42
    7.28 ,  15.26 ,  29.54 ,   0.   ,   0.   ,   0.   ,&               !43
    7.37 ,  16.76 ,  28.47 ,   0.   ,   0.   ,   0.   ,&               !44
    7.46 ,  18.08 ,  31.06 ,   0.   ,   0.   ,   0.   ,&               !45
    8.34 ,  19.43 ,  32.93 ,   0.   ,   0.   ,   0.   ,&               !46
    7.576,  21.49 ,  34.83 ,   0.   ,   0.   ,   0.   ,&               !47
    8.993,  16.908,  34.48 ,   0.   ,   0.   ,   0.   ,&               !48
    5.786,  18.869,  28.03 ,  54.0  ,   0.   ,   0.   ,&               !49
    7.344,  14.632,  30.502,  40.734,  72.28 ,   0.   ,&               !50
    8.641,  16.53 ,  25.3  ,  44.2  ,  56.0  , 108.0  ,&               !51
    9.009,  18.6  ,  27.96 ,  37.41 ,  58.75 ,  70.7  ,&               !52
   10.451,  19.131,  33.0  ,   0.   ,   0.   ,   0.   ,&               !53
   12.130,  21.21 ,  32.1  ,   0.   ,   0.   ,   0.   ,&               !54
    3.894,  25.1  ,  35.0  ,   0.   ,   0.   ,   0.   ,&               !55
    5.212,  10.004,  37.000,   0.   ,   0.   ,   0.   ,&               !56
    5.577,  11.06 ,  19.177,  49.954,   0.   ,   0.   /                !57
  DATA (PotIon(i),i=343,456)    /&
   5.47 ,  10.85 ,  20.197,  36.758,   0.   ,   0.   ,&               !58
   5.42 ,  10.55 ,  21.624,  38.981,   0.   ,   0.   ,&               !59
   5.49 ,  10.72 ,  22.14 ,  40.42 ,   0.   ,   0.   ,&               !60
   5.55 ,  10.90 ,  22.42 ,  41.09 ,   0.   ,   0.   ,&               !61
   5.63 ,  11.07 ,  23.45 ,  41.47 ,   0.   ,   0.   ,&               !62
   5.67 ,  11.25 ,  24.71 ,  42.65 ,   0.   ,   0.   ,&               !63
   6.14 ,  12.1  ,  20.38 ,  44.03 ,   0.   ,   0.   ,&               !64
   5.85 ,  11.52 ,  21.98 ,  39.84 ,   0.   ,   0.   ,&               !65
   5.93 ,  11.67 ,  22.83 ,  41.56 ,   0.   ,   0.   ,&               !66
   6.02 ,  11.80 ,  22.84 ,  42.51 ,   0.   ,   0.   ,&               !67
   6.10 ,  11.93 ,  22.74 ,  42.66 ,   0.   ,   0.   ,&               !68
   6.18 ,  12.05 ,  23.68 ,  42.69 ,   0.   ,   0.   ,&               !69
   6.254,  12.17 ,  25.03 ,  43.74 ,   0.   ,   0.   ,&               !70
   5.426,  13.9  ,  20.960,  45.193,   0.   ,   0.   ,&               !71
   7.0  ,  14.9  ,  23.3  ,  33.319,   0.   ,   0.   ,&               !72
   7.89 ,  16.200,  24.0  ,   0.   ,   0.   ,   0.   ,&               !73
   7.98 ,  17.70 ,  25.0  ,   0.   ,   0.   ,   0.   ,&               !74
   7.88 ,  16.6  ,  26.0  ,   0.   ,   0.   ,   0.   ,&               !75
   8.7  ,  17.0  ,  27.0  ,   0.   ,   0.   ,   0.   /                !76
  DATA (PotIon(i),i=457,594)    /&
   9.1  ,  20.0  ,  28.0  ,   0.   ,   0.   ,   0.   ,&               !77
   9.0  ,  18.563,  29.0  ,   0.   ,   0.   ,   0.   ,&               !78
   9.225,  20.5  ,  30.0  ,   0.   ,   0.   ,   0.   ,&               !79
  10.437,  18.756,  34.2  ,   0.   ,   0.   ,   0.   ,&               !80
   6.108,  20.428,  29.83 ,   0.   ,   0.   ,   0.   ,&               !81
   7.416,  15.032,  31.937,  42.32 ,  68.8  ,   0.   ,&               !82
   7.289,  16.69 ,  25.56 ,  45.3  ,  56.0  ,  88.3  ,&               !83
   8.42 ,  19.0  ,  27.0  ,   0.   ,   0.   ,   0.   ,&               !84
   9.3  ,  20.0  ,  30.   ,   0.   ,   0.   ,   0.   ,&               !85
  10.748,  20.0  ,  30.0  ,   0.   ,   0.   ,   0.   ,&               !86
   4.0  ,  22.0  ,  33.0  ,   0.   ,   0.   ,   0.   ,&               !87
   5.279,  10.147,  34.0  ,   0.   ,   0.   ,   0.   ,&               !88
   6.9  ,  12.1  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !89
   6.0  ,  11.5  ,  20.0  ,  28.8  ,   0.   ,   0.   ,&               !90
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !91
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !92
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !93
   5.800,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !94
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !95
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !96
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !97
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   ,&               !98
   6.0  ,  12.0  ,  20.0  ,   0.   ,   0.   ,   0.   /                !99
   save

 wl=w/10.d0
 elow=e*8066.d0
 epup=1.d7/wl+elow

   icode=code+0.5
  dde=(code-FLOAT(icode))*100.+1.01
  Zeff=dde
  Z2=Zeff*Zeff
  w2=wl*wl
  Elow=Elow*2.997925d10
  n1=code
  n2=dde+0.5
  ne=(n1-1)*6+n2
   if(gr == 0.d0)gr=log10(2.223d13/w2)
   if(gs == 0.d0)then
!    eff=25.
!    dEpup=PotIon(ne)-Epup*1.23981d-4
!    if(dEpup > 0.)eff=13.595d0*Z2/dEpup
!    E2=eff*eff
!    gs=1.d-8*E2*SQRT(eff)
! Same as in Piskunov's Synth code
     enu4=(Z2*13.598/(Potion(ne)-Epup*1.23981d-4))**2
     if(n2 == 1)then
       gqst=2.26d-7*enu4
     else
       gqst=5.42d-7*enu4/(n2+1)**2
     endif
     gs2=log10(gqst)
!    gs=log10(gs)
!    write(56,*)code,gs,gs2
    gs=gs2
   endif
   if(gw == 0.d0)then
    eff=25.
    dEpup=PotIon(ne)-Epup*1.23981d-4
    if(dEpup > 0.)eff=13.595d0*Z2/dEpup
    E2=eff*eff
    rsq=2.5d0*E2/Z2
    Nseq=code-E+1
    gw=4.5d-9*rsq**.4d0
    gw=log10(gw)
   endif

 end  
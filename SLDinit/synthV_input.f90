subroutine SynthV_input ! read input data and open files

use synth_data
use chemical_elements
use atmosphere_model
use mol_dat
use dflib
character(5) new_VALD
character(80) file_name,string,dbase_name

 open(2,file='SynthV.config',status='old',iostat=ios) ! Configuration file, contained all need parameters for calculation
 if(ios /= 0)stop 'Configuration file SynthV.config does not exist' 

! Lines lists 
 read(2,'((a))')file_name ! Name of file with atomic lines list
 k=index(file_name,'!')  ! Skip comments
 if(k == 0)k=81
 k=k-1
 m=index(file_name(1:k),'lns') ! Atomic lines in binary ".lns" format
 if(m /= 0)then
  open(3,file=file_name(1:k),form='binary',access='direct',recl=18,status='old',iostat=ios)
  if(ios == 0)atomic_lines=1
 else
  open(3,file=file_name(1:k),status='old',iostat=ios)
  if(ios == 0)then
   read(3,'((a))')string
   m=index(string,'Wavelength region')
   if(m /= 0)then
    atomic_lines=2  ! Atomic lines in ASCII VALD format
    read(string,*)d,d,nVALD ! Check if in lines list are molecular lines
	read(3,*)
	read(3,'((a))')string
	if_ev=1
	k=index(string,'(cm^')
    if(k /= 0)if_ev=0   ! Elow in kaisers
 	do i=1,nVALD
	 new_VALD=' '
	 read(3,*)new_VALD
	 if(molecular_lines == 2)cycle
     do l=1,Nmol
	  k=index(new_VALD,trim(Mol_id(l)))
	  if(k /= 0)then
       molecular_lines = 2; exit
      endif
     enddo
    enddo
   endif
   rewind (3)
  endif
 endif

 read(2,'((a))')file_name ! Name of file with molecular lines list
 k=index(file_name,'!')  ! Skip comments
 if(k == 0)k=81
 k=k-1
  open(13,file=file_name(1:k),form='binary',access='direct',recl=16,status='old',iostat=ios) ! molecular lines list in binary format
  if(ios == 0)molecular_lines=1
  
! Stellar atmosphere model
 file_name(1:80)=' '; model_type=0; model_abundance=abund
 read(2,'((a))')file_name  ! Name of file with model atmosphere
 k=index(file_name,'!')  ! Skip comments
 if(k == 0)k=81
 k=k-1
 open(1,file=file_name(1:k),status='old',iostat=ios)
 if(ios /= 0)stop 'Model atmosphere file does not exist'
 read(1,'((a))')string
 l=index(string,'LTE')
 if(l /= 0)then
  model_type=1 ! Kurucz's format atmosphere model
 else
  read(1,'((a))')string
  l=index(string,'T EFF')
  if(l /= 0)model_type=2 ! Piskunov's .krz format
  l=index(string,'Teff [K]')
  if(l /= 0)model_type=3 ! MARCS format
 endif
 if(model_type == 0)stop 'Unknown model atmosphere format'
 rewind (1)
 select case (model_type)  ! Load model atmosphere
  case (1)
   call load_modelK ! Kurucz's format
  case (2)
   call load_modelP !Piskunov's format
  case (3)
   call load_modelM ! MARCS format
 end select     

 read(2,*,iostat=ios)wavelength_first,wavelength_last,wavelength_step ! Describe the wavelengths region
 if(ios /= 0.or.wavelength_last < wavelength_first)stop 'Wrong wavelengths region record'
 if(atomic_lines == 1)call first_record(3,first_record_atomic,wavelength_first-2.5d0) ! first record in atomic lines list
 if(molecular_lines == 1)call first_record(13,first_record_molecular,wavelength_first-2.5d0) ! first record in molecular lines list
 number_of_wavelengths=(wavelength_last-wavelength_first)/wavelength_step+1

 read(2,*)V_micro ! Microturbulent velocity (km/sec)
 Vturb=V_micro    ! Fixed with depths

 file_name(1:80)=' '
 read(2,'((a))')file_name ! Name of output file
 k=index(file_name,'!')  ! Skip comments 
 if(k == 0)k=81
 k=k-1
 dbase_name=file_name(1:k)
 output_format='my '; l=index(file_name(1:k),'prf'); if(l /= 0)output_format='prf'; l=index(file_name(1:k),'out'); if(l /= 0)output_format='out'
  open(16,file=file_name(1:k)) ! ASCII output
  if(output_format == 'my ')then 
   k=index(file_name(1:k),'.',BACK=.TRUE.); file_name=file_name(1:k)//'bsn' ! Binary output file
   open(26,file=file_name, form='binary',access='direct',recl=144)
   write(26,rec=1)wavelength_first,wavelength_last,wavelength_step,number_of_wavelengths,Teff,glog,V_micro
   file_name=file_name(1:k)//'lin' ! Lines identification file
   open(27,file=file_name)
  endif

 file_name(1:80)=' '
 read(2,'((a))')file_name ! File name with individual abundancies
 k=index(file_name,'!')  ! Skip comments
 if(k /= 0)file_name(k:80)=' '
 open(100,file=file_name,status='old',iostat=ios)
 if(ios == 0)then
  do j=1,10    ! Read individual abundancies from external file
   read(100,*)
   ns=(j-1)*10+1
   nf=ns+9
   if(j == 10)nf=ns+8
    read(100,*)(abund(k),k=ns,nf)
  enddo
 else
  l=index(file_name,'solar');l1=index(file_name,'Solar');l2=index(file_name,'SOLAR') ! use build in solar abundance
  if(atomic_lines == 2)then
   call check_VALD_abundance
   rewind(3)
  endif
  if(l+l1+l2 == 0)abund=model_abundance
 endif
 close (100)

 file_name(1:80)=' '
 read(2,'((a))')file_name ! File name with data of abundancies stratification
 k=index(file_name,'!')  ! Skip comments
 if(k /= 0)file_name(k:80)=' '
 open(100,file=file_name,status='old',iostat=ios)
 if(ios == 0)call read_stratification
  close (100)

 file_name(1:80)=' '
 read(2,'((a))')file_name ! File name with Vturb vs. depths
 k=index(file_name,'!')  ! Skip comments
  if(k /= 0)file_name(k:80)=' '
 open(100,file=file_name,status='old',iostat=ios)
 if(ios == 0)then
  do j=1,Nrhox
   read(100,*)Vturb(j)
  enddo
 endif
 close (100)
 
 file_name(1:80)=' '; Vrad=0.d0; Vrad_switch=0
 read(2,'((a))')file_name ! File name with Vrad vs. depths
 k=index(file_name,'!')  ! Skip comments
 if(k /= 0)file_name(k:80)=' '
 open(100,file=file_name,status='old',iostat=ios)
 if(ios == 0)then
 read(100,*,iostat=ios)Vrad(1)
 if(ios == 0)then
  Vrad_switch=1
    do j=2,Nrhox
     read(100,*,iostat=ios)Vrad(j)
	 if(ios /= 0)then
	  Vrad(2:Nrhox)=Vrad(1); Vrad_switch=2; exit ! Fixed Vrad
     endif
    enddo
  else
   stop 'Wrong data in Vrad vs depths file'
  endif
 endif
 close (100) 

 length_record=8*Nrhox
 open(200,file='SynthV.tmp',form='binary',access='direct',recl=length_record) ! Temporary file for lines absorbtion

 if(atomic_lines + molecular_lines == 0)write(*,*)'Beware: both atomic and molecular lines has been not included!'

 read(2,'((a))',iostat=ios)string ! Calculate database of individual lines profiles
 if(ios == 0)then
  k=index(string,'profiles')
  if(k /= 0) if_database=1
  k=index(dbase_name,'.',BACK=.TRUE.); file_name=file_name(1:k)//'bsn'
  open(133,file=dbase_name(1:k)//'dtbs',form='binary') ! database with lines profile
  open(134,file=dbase_name(1:k)//'dir',form='binary') ! information about spectral lines in database
  proffilename = dbase_name(1:k)//'adp'
  proffilename_inf = dbase_name(1:k)//'dir'
  proffilename_bin = dbase_name(1:k)//'dbs'

  itest=0
  read(2,'((a))',iostat=ios) string
  if(ios == 0)then
 	k=index(string,'test')
	if(k /= 0) itest=1
  endif
 endif


end

subroutine calculate_tau5000 ! Calculate continuum opacity and optical depths for 5000 A (need for effective depths calculations)
use radiation_field
use synth_data
use atmosphere_model

  freq=2.997925d18/5000.d0
  freqlg=log(freq)
  do k=1,Nrhox
    EhvkT(K)=exp(-freq*hkT(k))
    stim(k)=1.d0-EhvkT(K)
    Bnu(k)=1.47439d-27*freq**3*EhvkT(k)/stim(k)*1.d-20
  enddo

  call KAPP
  Aline=0.d0
  call tau_nu

  tau5000=taunu
  kappa5000=acont+sigmac
  logtau5000=log10(tau5000)

end

subroutine read_stratification ! Read abundance stratification in same as in LLmodels format
use  synth_data
use atmosphere_model
use chemical_elements
character(300) string

 string(1:300)=' '
 read(100,'((a))')string
 number_stratified=0
  do i=1,99
   k=index(string, Elem_id(i))
   if(k == 0)cycle
   number_stratified=number_stratified+1
   atomic_numbers_stratified(number_stratified)=i
  enddo

  do j=1,Nrhox
   read(100,*)dummy,(xab_Strat(k,j),k=1,number_stratified)
   do k=1,number_stratified
    xab_Strat(k,j)=10**(xab_Strat(k,j))
   enddo
  enddo

end

subroutine check_VALD_abundance ! Read (if given) abundancies from VALD lines list file
use atmosphere_model
use chemical_elements
character(10) elem_string(99)

 read(3,*)d,d,Nlines; read(3,*);read(3,*)
 do i=1,Nlines
  read(3,*)
 enddo
 read(3,*)
 read(3,*,iostat=ios)elem_string(1)
 if( ios /= 0)return
 k=index(elem_string(1),'H :')
 if(k == 0)return
 backspace (3)
 read(3,*)elem_string
 do i=1,99
  read(elem_string(i)(4:10),*)model_abundance(i)
 enddo
 rewind (3)
end  


 
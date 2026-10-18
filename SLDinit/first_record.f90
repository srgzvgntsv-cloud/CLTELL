subroutine first_record(lunit,record_number,wavelength) ! Determine first line record number in the atomic or molecular lines list luinit
use DFPORT
use DFLIB
use synth_data, only : atomic_lines,molecular_lines
integer(4) stat_data(12),record_number
real(8) wavelength,wl1,wl2,wlmax
character(4)id_ml

 length_record=18
 if(lunit == 13)length_record=16 ! molecular lines list
 istat=FSTAT(lunit,stat_data)
 max_records=stat_data(8)/length_record ! Number of records in file

! Checking possibility to use internal lines lists for calculations in given wavelengths region 
 if(lunit == 13)then
  read(lunit,rec=max_records)id_ml,iw
  wlmax=iw/10000.d0
  if(wavelength > wlmax)then
   write(*,*)'Beware: molecular lines list is shorter then spectrum wavelengths'
   molecular_lines=0
   return
  endif
 endif  

 i=1; j=max_records+1 ! Quick search algorithm
 do
  k=(i+j)/2
  if(k == 0)k=1;
  if(k > max_records-1)then
   j=k-1; i=1
  else 
   if(lunit == 3)then ! Atomic
    read(lunit,rec=k)iw
	wl1=iw/1000.d0
    read(lunit,rec=k+1)iw
    wl2=iw/1000.d0
   else
    read(lunit,rec=k)id_ml,iw ! molecular
    wl1=iw/10000.d0
    read(lunit,rec=k+1)id_ml,iw
    wl2=iw/10000.d0
   endif
   if(wavelength >= wl1.and.wavelength <= wl2)then
    record_number=k
    exit
   endif
   if(wavelength < wl1)j=k
   if(wavelength > wl1)i=k
  endif
 enddo

end
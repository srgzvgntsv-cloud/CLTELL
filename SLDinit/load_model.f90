subroutine load_model ! Load model atmosphere file 

use atmosphere_model
use chemical_elements
character(80) string
real(8) pg,prad,pturb,pe 
real(8) abund_M(92)

entry load_modelK    ! Kurucz's (or LbL) format
 do
  read(1,'((a))')string
  k=index(string,'BEGIN')
  if(k /= 0)exit ! End of model atmosphere

  k=index(string,'TEFF')
   if(k /= 0)read(string(k+4:k+11),*)Teff
  k=index(string,'GRAVIT')
  if(k /= 0)read(string(k+7:k+12),*)glog
  k=index(string,'TITLE')
  if(k /= 0)read(string(k+5:80),'(a74)')Title
  k=index(string,'SCALE'); k1=index(string,'CHANGE')
  if(k /= 0.and.k1 /= 0)then
   backspace 1
   read(1,'(16x,f9.5,17x,2(i2,f8.5)/(17x,6(i3,f7.2)))')xscale,(k,model_abundance(i),i=1,99)
  scale_log=log10(xscale)
  else
   if(k /= 0)then
    backspace 1
    read(1,'(16x,f9.5)')xscale
    scale_log=log10(xscale)
   endif
  endif
  k=index(string,'READ DECK')
  if(k /= 0)then
   l=index(string,'READ DECK6')
   n=k+9; if(l /= 0)n=l+10
   do i=n,80
    if(string(i:i) /= ' ')exit
   enddo
   n1=i
   do i=n1,80
    if(string(i:i) == ' ')exit
   enddo
   n2=i
   read(string(n1:n2),*)Nrhox  ! Number of depths
   call arrays_allocation(Nrhox) ! Allocate arrays dependend from number of depths
     do j=1,Nrhox
      read(1,*)Rhox(j),T(j),P(j),Xne(j),abRoss(j)
     enddo
   endif
 enddo
 Grav=10**(Glog)
 do i=3,99
  model_abundance(i)=model_abundance(i)+scale_log
 enddo
 return

entry load_modelP ! Piskunov's .krz format

 read(1,'((a))',iostat=ios) string
 if( ios/=0 ) stop 'Error: Corrupted KRZ model atmosphere file'
!    Title = string(k+5:64)
    read(1,'((a))') string
	k1=index(string,'T EFF=');k2=index(string,'GRAV=');	k3=index(string,'MODEL')
	read(string(k1+6:k2-1),*)Teff; read(string(k2+5:k3-1),*)glog ! Beware - check if model is really motype == 1 or 3
!  read(1,'(2(6x,F6.0),12x,i1,8x,f8.1)') Teff,glog,motype,wl_std 
  read(1,*) ! Skip ifop data
  read(1,*) model_abundance,Nrhox
  call arrays_allocation(Nrhox) ! Allocate arrays dependent from number of depths
  do j=1,nrhox
    read(1,*,iostat=ios)Rhox(j),T(j),xne(j),xnatom(j),Rho(j)
    P(j)=(xnatom(j)+xne(j))*1.38054d-16*T(J) 
  enddo
 Grav=10**(Glog)
 return

entry load_modelM ! MARCS models format

do
 string=' '
 read(1,'((a))')string
 k=index(string,'Teff'); if(k /= 0)read(string,*)Teff
 k=index(string,'gravity')
  if(k /= 0)then
   read(string,*)Grav; Glog=log10(Grav);  write(Title,"('MARCS model Teff=',f7.0,' logg=',f6.2)")Teff,Glog
  endif
 k=index(string,'Logarithmic')
  if(k /= 0)then							  
   read(1,*)abund_m
    do i=1,92
     if(abund_m(i) < -90)then
      model_abundance(i)=abund_m(i)
     else                                           ! Transform abundancies from 12.0 log scale to 1.0 scale
      model_abundance(i)=10**(abund_m(i)-abund_m(1))*abund(1)
      if(i > 2)model_abundance(i)=log10(model_abundance(i))
     endif
    enddo
    model_abundance(93:99)=abund(93:99)
  endif
 k=index(string,'Number of depth points')
 if(k /= 0)then
  read(string,*)Nrhox
  call arrays_allocation(Nrhox) ! Allocate arrays dependent from number of depths
 endif
 k=index(string,' k lgTauR')
 if(k /= 0)then
  do j=1,Nrhox
   read(1,*)k,tau,taus,z,T(j),pe,pg,prad,pturb
   P(j)=pg+prad+pturb
   Xne(j)=pe/ (1.38054d-16 * T(j))
  enddo
  read(1,*)
  do j=1,Nrhox
   read(1,*) k,taudum,abRoss(j),ro,emu,v,frconv,Rhox(j)
  enddo
  exit
 endif
enddo    
return

end
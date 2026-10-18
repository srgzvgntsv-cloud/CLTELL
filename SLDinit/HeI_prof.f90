subroutine He_prof(i_He,j,dw,pl) ! Calculate profile for one of 20 HeI lines:
 use atmosphere_model
 use lines_data
 use HeTab
 implicit real(8) (a-h,o-z)
 integer(2), save :: first_time =0
 integer(2) i_He
 integer(4) j
 real(8) wtmp(4)
  if(first_time == 0)then  ! Initialize Griem's constants
   first_time=1;    call He_Griem

  endif
!!!!!!!!!!!!!!!!!!
! else
 if(i_He <= 20)then
  omeg=omega(i_He,j)
  alph=alp(i_He,j)
  dgrm=dgr(i_He,j)
  TT=T(j)
  xe=xNe(j)
  wl=central_wavelength
  dump=Prof(omeg,alph,dgrm,TT,xe,dw,d,wl)
   ww=dw+d
   admp=adump(j)+dump/10.d0*Dopdump(j)
   ww=abs(ww)/10.d0
   xDop=ww*Dopdump(j)
   vv=Voigt(xDop,admp)
endif
!!!!!!!!!!!!!!!!!!

 vvtemp=0.d0
 if(i_He == 9.and.xNe(j) >= 1.d13.and.xNe(j) <=2.138d16.and.T(j) <= 40000.)then
! 4471
   vvtmp=He4471(T(j),xNe(j),dw)
   if(vvtmp > 0.0) vv=vvtmp
!!!!!!!!!!!!!!!!!!
 else if(i_He == 4.and.xNe(j) > 1.d14.and.xNe(j) <=6.d17.and.T(j) <= 40000.0.and.T(j) >= 10000.)then
! 4026
   vvtmp=He_pr(T(j),xNe(j),dw,1)
   if(vvtmp > 0.0) vv=vvtmp

 else if(i_He == 22.and.xNe(j) > 1.d14.and.xNe(j) <=6.d17.and.T(j) <= 40000.0.and.T(j) >= 10000.)then
! 4109
   vvtmp=He_pr(T(j),xNe(j),dw,2)
   if(vvtmp > 0.0) vv=vvtmp

 else if(i_He == 23.and.xNe(j) > 1.d14.and.xNe(j) <=6.d17.and.T(j) <= 40000.0.and.T(j) >= 10000.)then
! 4143
   vvtmp=He_pr(T(j),xNe(j),dw,3)
   if(vvtmp > 0.0) vv=vvtmp

 else if(i_He == 6.and.xNe(j) > 1.d14.and.xNe(j) <=6.d17.and.T(j) <= 40000.0.and.T(j) >= 10000.)then
! 4387
   vvtmp=He_pr(T(j),xNe(j),dw,4)
   if(vvtmp > 0.0) vv=vvtmp
 else if(i_He == 10.and.xNe(j) > 1.d14.and.xNe(j) <=4.677E16.and.T(j) <= 40000.and.dw > -2.1d0)then
! 4922
   vvtmp=He_pr(T(j),xNe(j),dw,5)
   if(vvtmp > 0.0) vv=vvtmp

!!!!!!!!!!!!!!!!!!
 else if (i_He>20.and.vvtemp <= 0.d0) then
   ii=i_He-21
   d=0.0

     do k=1,4
      wtmp(k)=wwe(k,ii)
     enddo
   id=MAP1(tgr,wtmp,4,T(J),ade,1)
   dump=ade*xNe(j)*1.d-16
     do k=1,4
      wtmp(k)=wwi(k,ii)
     enddo
   id=MAP1(tgr,wtmp,4,T(J),ade,1)
   dump=dump+ade*xDens(j,1,2)*1.d-16
   ww=dw+d
   admp=adump(j)+dump/10.d0*Dopdump(j)
   ww=abs(ww)/10.d0
   xDop=ww*Dopdump(j)
   vv=Voigt(xDop,admp)
!!!!!!!!!!!!!!!!!!
 endif

 pl = al_int(j)*vv

end



subroutine He_Griem ! Interpolate HeI constants from Griem textbook
 use atmosphere_model
 use HeTab
 implicit real(8) (a-h,o-z)
 real(8) bwid(4),alpp(4),sh(4),shft(1)
   do i=1,20
     do k=1,4
      bwid(k)=bw(k,i)
      alpp(k)=alphg(k,i)
      sh(k)=shift(k,i)
     enddo
      do  j=1,Nrhox
       id=MAP1(tgr,bwid,4,T(J),omega(i,j),1)
       id=MAP1(tgr,alpp,4,T(J),alp(i,j),1)
       id=MAP1(tgr,sh,4,T(J),shft,1)
       dgr(i,j)=shft(1)/omega(i,j)
       omega(i,j)=omega(i,j)*xNe(j)*1.d-16
       alp(i,J)=alp(i,J)*SQRT(SQRT(xNe(j)))*1.d-4
      enddo
    enddo
end



real(8) function Prof(omeg,alph,di,T,xe,dw,d,w) ! Stark profile from UT
 use HeTab
 implicit real(8) (a-h,o-z)
 real(8) wt(18),dt(18),k
  x=dw/omeg-di
  ro0=(2.3873241d-1/xe)**0.3333333333333d0
  v=1.432d6*SQRT(T/1.d4)
  we=1.88365d19/(w*w)*omeg
  rd=6.9d0*SQRT(T/xe)
  sig=we*ro0/v
  as=alph**0.8888888888889d0/sig**0.33333333333333333333d0
  k=1.03d0*sig**1.33333333333d0*alph**0.44444444444444d0*x
  i=2
  if(k > 0.)i=1
  sgn=1.d0
  if(alph < 0.)SGN=-1.d0
  ak=ABS(k)
   if(ak > 5.)then
     tk=SQRT(SQRT(ak))
    if(i == 1)then
     wr=0.d0
     dr=3.31d0*tk
    endif
    if(i == 2)then
     wr=2.34d0*tk
     dr=wr
    endif
   else
    do j=1,18
     wt(j)=wi(i,j)
     dt(j)=dik(i,j)
    enddo
    id=MAP1(kt,wt,18,ak,wr,1)
    id=MAP1(kt,dt,18,ak,dr,1)
   endif
    d=sgn*(as*dr-3.d0*alph**1.33333333333d0*ro0/rd)
    d=-omeg*(di+d)
     Prof=(1.d0+as*wr)*omeg
end



real(8) function He4471(Tj,xj,dw) ! Interpolate Stark profile for HeI 4471A
 use HeTab
 implicit real(8) (a-h,o-z)
  do j=1,7
   if(ne(j) >= xj)then
    nne=j
    exit
   endif
  enddo
  do j=1,4
   if(Tgr(j) >= Tj)then
    nT=j
    exit
   endif
  enddo
  nn1=nne-1
  nT1=nT-1
  n11=n0(nn1)
  nd1=nd(nn1)+n11-1
  n12=n0(nne)
  nd2=nd(nne)+n12-1
  m2=0
   do j=n11,nd1
    if(a(j) >  dw)then
     m2=j
     exit
    endif
   enddo
!!!!!!!!!!!!!!!!!!
!  if(m2 == 0.or.m2 == 1) then
!    He4471=0.0
!    goto 555
!  endif
  if(m2 == 0)m2=nd1
  m1=m2-1
   if(m1 <  n11)then
    m1=m2
    m2=m2+1
   endif
   l2=0
    do j=n12,nd2
     if(a(j) > dw)then
      l2=j
      exit
     endif
    enddo
!!!!!!!!!!!!!!!!!!
!  if(l2 == 0.or.l2 == 1) then
!    He4471=0.0
!    goto 555
! endif
  if(l2 == 0)l2=nd2
   l1=l2-1
    if(l1 <  N12)then
     l1=l2
     l2=l2+1
    endif
   dt=Tgr(nt)-Tgr(NT1)
   dne=ne(nne)-ne(nn1)
   dt=(Tj-Tgr(nt1))/dt
   dne=(xj-ne(nn1))/dne
   dw1=(dw-a(m1))/(a(m2)-a(m1))
   dw2=(dw-a(l1))/(a(l2)-a(l1))
   ds11=Sd(nt1,m1)+dt*(Sd(nt,m1)-Sd(nt1,m1))
   ds12=Sd(nt1,m2)+dt*(Sd(nt,m2)-Sd(nt1,m2))
   ds21=Sd(nt1,l1)+dt*(Sd(nt,l1)-Sd(nt1,l1))
   ds22=Sd(nt1,l2)+dt*(Sd(nt,l2)-Sd(nt1,l2))
    ds1=ds11+dw1*(ds12-ds11)
    ds2=ds21+dw2*(ds22-ds21)
      He4471=ds1+dne*(ds2-ds1)*0.56418958d0
      ddop=9.6167475d-2*SQRT(Tj/1.d4)
      He4471=He4471*ddop*1.7724538d0
555 continue
end
real(8) function He_pr(Tj,xj,dw,l_He) 
!
! Interpolate Stark profile for HeI 4922A
! Profiles for HeI 4922A line from 12.A&A.542.A75 with Doppler conv:
!
! 4026, 4109, 4143, 4387 A
! Tabulated HeI line profiles after Beauchamp et al (ApJS, 108, 559, 1997)
! Copied from SPECTRUM v2.72 code (see http://www.phys.appstate.edu/spectrum/spectrum.html)
!


 use HeTab
 implicit real(8) (a-h,o-z)
  do j=1,9
   if(ne49(j,l_He) >= xj)then
    nne=j
    exit
   endif
  enddo
  do j=1,4
   if(Tgr(j) >= Tj)then
    nT=j
    exit
   endif
  enddo
  nn1=nne-1
  nT1=nT-1
  n11=n049(nn1,l_He)
  nd1=nd49(nn1,l_He)+n11-1
  n12=n049(nne,l_He)
  nd2=nd49(nne,l_He)+n12-1
  m2=0

   do j=n11,nd1
    if(a49(j,l_He) >  dw)then
     m2=j
     exit
    endif
   enddo
!!!!!!!!!!!!!!!!!!
  if(m2 == 0.or.m2 == 1) then
    He_pr=0.0
    goto 555
  endif
  m1=m2-1
   if(m1 <  n11)then
    m1=m2
    m2=m2+1
   endif
   l2=0
    do j=n12,nd2
     if(a49(j,l_He) > dw)then
      l2=j
      exit
     endif
    enddo
!!!!!!!!!!!!!!!!!!
  if(l2 == 0.or.l2 == 1) then
    He_pr=0.0
    goto 555
  endif
   l1=l2-1
    if(l1 <  N12)then
     l1=l2
     l2=l2+1
    endif
   dt=Tgr(nt)-Tgr(NT1)
   dne=ne49(nne,l_He)-ne49(nn1,l_He)
   dt=(Tj-Tgr(nt1))/dt
   dne=(xj-ne49(nn1,l_He))/dne
   dw1=(dw-a49(m1,l_He))/(a49(m2,l_He)-a49(m1,l_He))
   dw2=(dw-a49(l1,l_He))/(a49(l2,l_He)-a49(l1,l_He))
   ds11=Sd49(nt1,m1,l_He)+dt*(Sd49(nt,m1,l_He)-Sd49(nt1,m1,l_He))
   ds12=Sd49(nt1,m2,l_He)+dt*(Sd49(nt,m2,l_He)-Sd49(nt1,m2,l_He))
   ds21=Sd49(nt1,l1,l_He)+dt*(Sd49(nt,l1,l_He)-Sd49(nt1,l1,l_He))
   ds22=Sd49(nt1,l2,l_He)+dt*(Sd49(nt,l2,l_He)-Sd49(nt1,l2,l_He))
    ds1=ds11+dw1*(ds12-ds11)
    ds2=ds21+dw2*(ds22-ds21)
      He_pr=ds1+dne*(ds2-ds1)*0.56418958d0
      ddop=9.6167475d-2*SQRT(Tj/1.d4)
      He_pr=He_pr*ddop*1.7724538d0
555 continue
end

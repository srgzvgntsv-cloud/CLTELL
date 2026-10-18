SUBROUTINE JOSH ! Solve transfer equition by Harvard method

! Calculate of opacity depth's
 call Tau_nu
! Solve transfer equition :
 call transfer
! Calculate radiation field:
 call radfield
end

subroutine Tau_nu ! Opacity depth
 use radiation_field
 use atmosphere_model
 implicit real(8) (a-h,o-z)
! Total opacity:
 Abtot=Acont+Aline+Sigmac; Alpha=Sigmac/Abtot
 CALL INTEG(RHOX,ABTOT,TAUNU,NRHOX,ABTOT(1)*RHOX(1))

end

subroutine Transfer ! Solve trahsfer equition
 use radiation_field
 use atmosphere_model
 use JoshTab
 IMPLICIT REAL(8) (A-H,O-Z)
 real(8) xsbar(51),xalpha(51),diag(51)
 maxJ=0
 if(Taunu(1) > xtau(Nxtau))MaxJ=1
  IF(maxJ /= 1)then
! Interpolation on the "standard" set of tau's:
    maxJ=MAP1(TAUNU,Bnu,NRHOX,XTAU,XSBAR,NXTAU)
    maxJ=MAP1(TAUNU,ALPHA,NRHOX,XTAU,XALPHA,NXTAU)
! Check interpolation:
   do l=1,Nxtau
    xalpha(l)=max(xalpha(l),0.d0)
    xsbar(l)=max(xsbar(l),1.d-38)
     if(xtau(l) < taunu(1))then
      xsbar(l)=max(bnu(1),1.d-38)
      xalpha(l)=max(alpha(1),0.d0)
     endif
    xs(l)=xsbar(l)
    diag(l)=1.d0-xalpha(l)*coefJ(l,l)
    xsbar(l)=(1.d0-xalpha(l))*xsbar(l)
   enddo
! Solution of Milne-Eddington equition: maximum number of iteration is nxtau
   DO L=1,NXTAU
    IFERR=0
    K=NXTAU+1
     DO KK=1,NXTAU
      K=K-1
      DELXS=0.d0
       DO  M=1,NXTAU
        DELXS=DELXS+COEFJ(K,M)*XS(M)
       enddo
      DELXS=(DELXS*XALPHA(K)+XSBAR(K)-XS(K))/DIAG(K)
      ERRORX=ABS(DELXS/XS(K))
      IF(ERRORX > .00001)IFERR=1
      XS(K)=MAX(XS(K)+DELXS,1.d-38)
     enddo
     IF(IFERR == 0)exit
    enddo
! Backward interpolation of source function:
    mdummy=MAP1(xtau,xs,nxtau,taunu,Snu,maxj)
  endif
 if(maxJ == Nrhox)return
! Approach to the large tau's:
 maxJ1=maxJ+1
 if(maxJ == 1)maxJ1=1
  do j=maxJ1,Nrhox
   Snu(j)=Bnu(j)
  enddo
 m=max(maxJ-1,1)
 nm1=Nrhox-m+1
 nmj=Nrhox-maxJ+1
! Solution of Milne-Eddington equition:
 do l=1,nxtau
  error=0.d0
  call DERIV(Taunu(m),Snu(m),Hnu(m),nm1)
   do j=m,Nrhox
    Hnu(j)=Hnu(j)/3.d0
   enddo
  call DERIV(taunu(maxJ),Hnu(maxJ),Jmins(maxj),nmj)
   do j=maxJ1,Nrhox
    Jnu(j)=Jmins(j)+Snu(j)
    Snew=(1.d0-alpha(j))*Bnu(j)+alpha(j)*Jnu(j)
    error=abs(Snew-Snu(j))/Snew+error
    Snu(j)=Snew
   enddo
  if(error < 0.00001)exit
 enddo
end

subroutine radfield ! Calculate Jnu,Hnu and specific intensities for fixed mu's
 use atmosphere_model
 use radiation_field
 use temporary_arrays
 use JoshTab
 
 implicit real(8) (A-H,O-Z)
 real(8) :: NEW,surft(7),EXTAU(51,7)
 data EXTAU /357*0.d0/

 if( maxJ ==1) then ! Specific intensity for very large tau (Jnu and Hnu calculated in Transfer)
  call PARCOE(Snu,Taunu,Aa,Bb,Cc,Nrhox)
  n1=Nrhox-1
   do j=1,Nrhox
    Ctwo(j)=Cc(j)*2.d0
    B2ct(j)=Bb(j)+Ctwo(j)*Taunu(j)
   enddo
   do j=1,n1
    B2ct1(j)=Bb(j)+Ctwo(j)*Taunu(j+1)
   enddo
   do mu=1,7
    ANGLE=mu_fix(mu)
    old=1.d0
    sum=0.d0
     do j=1,n1
      id=0
      tangle=Taunu(j+1)/ANGLE
      new=EXP(-tangle)
      d=tangle-Taunu(j)/ANGLE
       if(d > 0.03)then
      SUM=SUM+OLD*(SNU(J)+(B2CT(J)+CTWO(J)*ANGLE)*ANGLE)- &
          NEW*(SNU(J+1)+(B2CT1(J)+CTWO(J)*ANGLE)*ANGLE)
       IF(TANGLE < 300.)then
        old=new
        cycle
       else
        SURFt(mu)=SUM
        id=1
        exit
       endif
     else
      DDDDD=1.d0
      IF(D >0.001)DDDDD=((((D/9.d0+1.d0)*D/8.d0+1.d0)*D/7.d0+1.d0)*&
                            D/6.d0+1.d0)*D/5.d0+1.d0
      SUM=SUM+NEW*(SNU(J)+(SNU(J)+B2CT(J)*ANGLE+(SNU(J)+(B2CT(J)+&
          CTWO(J)*ANGLE)*ANGLE)*(DDDDD*D/4.d0+1.d0)*D/3.d0)*D/2.d0)*D
      OLD=NEW
     endif
    enddo
     if( id == 0)SURFt(mu)=SUM+OLD*(SNU(NRHOX)+(B2CT(NRHOX)+CTWO(NRHOX)&
                            *ANGLE)*ANGLE)
    Inu(mu)=surft(mu)
   enddo

  return
 endif

  do l=1,Nxtau  ! "Normal" opacity
   xjs(l)=-xs(l) ! Calculate Jnu, Hnu
    do  m=1,Nxtau
     XJS(L)=XJS(L)+COEFJ(L,M)*XS(M)
   enddo
   XH(L)=0.d0
    DO M=1,NXTAU
     if(l /= 1)then
      XH(L)=XH(L)+COEFH(L,M)*XS(M)
     else
      xh(1)=xh(1)+ch(m)*xs(m)
     endif
    enddo
   enddo
    MDUMMY=MAP1(XTAU,XJS,NXTAU,TAUNU,JMINS,MAXJ) ! Interpolate Jnu from fixed tau's
    MDUMMY=MAP1(XTAU,XH,NXTAU,TAUNU,HNU,MAXJ)    ! Interpolate Hnu from fixed tau's
   XK=0.d0
    DO  M=1,NXTAU
     XK=XK+CK(M)*XS(M)
    enddo
   KNU(1)=XK
    DO J=1,MAXJ
     JNU(J)=JMINS(J)+SNU(J)
    enddo
  
 
 CALL PARCOE(XS,XTAU,Aa,Bb,Cc,NXTAU)   ! Specific intensity 
   N1=NXTAU-1
    DO J=1,NXTAU
     CTWO(J)=Cc(J)*2.d0
     B2CT(J)=Bb(J)+CTWO(J)*XTAU(J)
    enddo
    DO J=1,N1
     B2CT1(J)=Bb(J)+CTWO(J)*XTAU(J+1)
    enddo
 
 IF(EXTAU(1,1) == 0.d0)THEN
    do  mu=1,7
     ANGLE=mu_fix(mu)
      DO J=1,NXTAU
       TANGLE=XTAU(J)/ANGLE
       IF(TANGLE < 300.)EXTAU(J,mu)=EXP(-TANGLE)
      enddo
    enddo
  ENDIF
  
   do mu=1,7
    Inu(mu)=0.d0
    ANGLE=mu_fix(mu)
    SURFt(mu)=0.d0
     DO J=1,N1
      IF(EXTAU(J,mu) == 0.d0)exit
      SURFt(mu)=SURFt(mu)+EXTAU(J,mu)*(XS(J)+(B2CT(J)+CTWO(J)*ANGLE)*ANGLE)-&
                EXTAU(J+1,mu)*(XS(J+1)+(B2CT1(J)+CTWO(J)*ANGLE)*ANGLE)
     enddo
      if(EXTAU(J,mu) /= 0.d0)then
       SURFt(mu)=SURFt(mu)+EXTAU(NXTAU,mu)*(XS(NXTAU)+(B2CT(NXTAU)+CTWO(NXTAU)*ANGLE)*ANGLE)
       Inu(mu)=surft(mu)
      endif
   enddo

end

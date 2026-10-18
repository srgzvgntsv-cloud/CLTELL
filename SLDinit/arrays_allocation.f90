subroutine arrays_allocation(n) ! Allocate all arrays with number of depths in dimensions
 use synth_data
 use atmosphere_model
 use radiation_field
 use lines_data 
 use temporary_arrays
 use mol_dat
 use HeTab
 use HlineTab

 allocate (T(n),Tkev(n),Tk(n),hkT(n),Tlog(n),abRoss(n),tauros(n),stat=j1)
 allocate (Rhox(n),xDens(n,99,6),xDens2(n,99,6),P(n),xNe(n),stat=j2)
 allocate (xNatom(n),rho(n),EhvkT(n),Stim(n),Bnu(n),stat=j3)
 allocate (Acont(n),Scont(n),Aline(n),Sline(n),sigmaC(n),stat=j4)
 allocate (Abtot(n),Alpha(n),AHline(n),Vturb(n),Vrad(n),stat=j5)
 allocate (Taunu(n),Snu(n),Hnu(n),Jnu(n),Jmins(n),Knu(n),stat=j6)
 allocate (al_int(n),Adump(n),Dopdump(n),stat=j7)
 allocate (xnmol(n,ndim),x_mol(n,nmol),omega(20,n),alp(20,n),dgr(20,n),stat=j9)
 j=j1+j2+j3+j4+j5+j6+j7+j9

 allocate (frac(n,6),PFPLUS(n),PFMIN(n),Told(n),Aa(n),Bb(n),Cc(n),stat=i1)
 allocate (xh(n),xjs(n),CTWO(n),B2CT(n),B2CT1(n),THETA(n),FFTHETA(n),stat=i2)
 allocate (BOLT(n,10),EXLIM(n),BOLTN(n,27),BOLTEX(n),FREET(n),w(n),dt(n),stat=i3)
 allocate (C1218(n),C1420(n),xHmin(n),C1240(n),C1444(n),bolt2(48,n),stat=i4)
 allocate (C1130(n),C1020(n),C1169(n),ahot(n),xnz(n,25),answer(n,6),stat=i5)
 allocate (E1(n,6),E2(n,6),xab_Strat(99,n),Tx(n),vt2(n),nt(n),H_cont(n),S_cont(n),J_cont(n),stat=i6)
 allocate (ProfIn(nHydSeries,55,nHydLines,n),xHh(nHydSeries,n),Mlast(n),stat=i7)
 allocate (aline_buffer(n),xnebuf(n),tau5000(n),kappa5000(n),logtau5000(n),stat=i8)
 i=i1+i2+i3+i4+i5+i6+i7+i8

 if(j+i /= 0) stop 'Memory allocation failed'
 return

 entry arrays_deallocation
  deallocate (T,Tkev,Tk,hkT,Tlog,abRoss,tauros,Rhox,xDens,xDens2,P,xNe,xNatom,rho,EhvkT,Stim,Bnu,&
  Acont,Scont,Aline,Sline,sigmaC,Abtot,Alpha,AHline,Vturb,Vrad,Taunu,Snu,Hnu,Jnu,Jmins,Knu,al_int,Adump,Dopdump,&
  xnmol,x_mol,omega,alp,dgr,frac,PFPLUS,PFMIN,Told,Aa,Bb,Cc,xh,xjs,CTWO,B2CT,B2CT1,THETA,FFTHETA,&
  BOLT,EXLIM,BOLTN,BOLTEX,FREET,w,dt,C1218,C1420,xHmin,C1240,C1444,bolt2,C1130,C1020,C1169,ahot,xnz,answer,&
  E1,E2,xab_Strat,Tx,vt2,nt,ProfIn,xHh,Mlast,xnebuf,tau5000,kappa5000,aline_buffer)

end

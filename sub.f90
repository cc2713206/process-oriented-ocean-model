MODULE sub    
USE param
CONTAINS

!******************
SUBROUTINE init
implicit none
INTEGER :: nb
real(8) :: di

ri(1)=r1; ri(2)=r2; ri(3)=r3; ri(4)=r4; ri(5)=r5 
ri(6)=r6; ri(7)=r7; ri(8)=r8; ri(9)=r9

do i =1, ni
  di=2*ri(i)
  if(di<=1.27) wi(i)=2.025*di**1.621
  if(di>1.27 .and. di<=7.) wi(i)=-0.103*di**2+4.069*di-2.024
  wi(i)=wi(i)/1000.
  if(ri(i)<=1./13.8**0.5)   nu(i)=1.+0.17*ri(i)*13.8**0.5
  if(ri(i)>1./13.8**0.5)  nu(i)=1.+0.55*ri(i)**(2./3.)*13.8**(1./3.)
  if(ri(i)<=1./2432.**0.5)   sh(i)=1.+0.17*ri(i)*2432.**0.5
  if(ri(i)>1./2432.**0.5)  sh(i)=1.+0.55*ri(i)**(2./3.)*2432.**(1./3.)
  re(i)=(1.5*ar)**(1./3.)*ri(i)/1000.
  uc2(i)=theta*(rho0-rhoi)*gv*2.*re(i)/rho0
  if(2.*re(i)<=0.001)then
    uc2(i)=(3.724*(2.*re(i)*100.)**0.256/100.)**2
  else
    uc2(i)=(7.656*(2.*re(i)*100.)**0.569/100.)**2  
  endif
  vi(i)=pi*ri(i)**3*ar*2.  
enddo

do i=1, ni-1
  dv(i)= vi(i+1)-vi(i)
enddo

tmix=smix*l1+l2+l3*(isdp+(isdp2-isdp)*0.)*dz

tpre=0.; tfre=0.

DO i = 0,nz+1
DO k = 0,nx+1
 p(i,k) = 0.0
 q(i,k) = 0.0
 dq(i,k) = 0.0
 u(i,k) = 0.0
 un(i,k) = 0.0
 ustar(i,k) = 0.0
 
 v(i,k) = 0.0
 vn(i,k) = 0.0
 vstar(i,k) = 0.0 
 
 w(i,k) = 0.0
 wn(i,k) = 0.0
 wstar(i,k) = 0.0
 wet(i,k) = .true.
 tfc(i,k)=0.
 us(i,k) = 0.   
 do j= 1, ni
   fc(j,i,k) = 0.    
 enddo
 az(i,k)=azmin
 taz(i,k)=azmin
 ek(i,k)=ekmin
 ep(i,k)=epmin 
 lzi=i*dz
 if(lzi<200.)then
   t(i,k)=(-1.92+1.7)/200.*lzi-1.7
 else
   t(i,k)=-1.92  
 endif   
 if(lzi<=90.)then
   s(i,k)=(34.7-34.22)/90.*lzi+34.22
 elseif(lzi>90. .and. lzi<=400.)then
   s(i,k)=(34.795-34.7)/(400.-90.)*(lzi-90.)+34.7
 else
   s(i,k)=(34.82-34.795)/(800.-400.)*(lzi-400.)+34.795
 endif  
 rho(i,k)=rho0*(1.+bs*(s(i,k)-s0)-bt*(t(i,k)-t0)) 
 rho(i,k)=rho(i,k)+tfc(i,k)*(rhoi-rho(i,k))
END DO
END DO

open ( 1, file = 'T0.dat' )
write(1,*) "variables=x,z,T"
write(1,*) "zone i=", nx, "j=", nz, "f=point"
DO i = 1, nz
  do k= 1, nx
    WRITE(1,'(f8.1,x,f6.1,x,f20.15,x,f20.15)') (k-1)*dx, -(i-1)*dz, t(i,k)
    enddo
  END DO  
close(1) 
open ( 1, file = 'S0.dat' )
write(1,*) "variables=x,z,S"
write(1,*) "zone i=", nx, "j=", nz, "f=point"
DO i = 1, nz
  do k= 1, nx
    WRITE(1,'(f8.1,x,f6.1,x,f20.15,x,f20.15)') (k-1)*dx, -(i-1)*dz, s(i,k)
    enddo
  END DO  
close(1) 

DO k = 0,nx+1
  usum(k) = 0.0
END DO

drdx = 1.0/(rho0*dx)
drdz = 1.0/(rho0*dz)

alpha=dt* 2.*om*sin(tta*pi/180.)

DO k = 0,nx+1
  nb=nz
  DO i = nb+1,nz+1
    wet(i,k) = .false.
  END DO
  if(k==0) then 
    do i=0, isdp
      wet(i,k)= .false.  
    enddo
  endif
END DO

DO k = 0,nx+1
DO i = 0,nz+1
  dry(i,k)=.NOT.wet(i,k)
END DO
END DO

DO k = 1,nx
DO i = 1,nz
  at(i,k) = dx/dz
  ab(i,k) = dx/dz
  ae(i,k) = dz/dx
  aw(i,k) = dz/dx
  IF(dry(i,k+1)) ae(i,k) = 0.0
  IF(dry(i,k-1)) aw(i,k) = 0.0
  IF(dry(i+1,k)) ab(i,k) = 0.0
  IF(dry(i-1,k)) at(i,k) = 0.0
  atot(i,k) = ab(i,k)+at(i,k)+ae(i,k)+aw(i,k)
END DO
END DO

ah=1.
aht=0.1
END SUBROUTINE init

!*******************
SUBROUTINE dyn
implicit none
real(8) :: pressx, pressz
real(8) :: perr, q1, q2, term1
INTEGER :: nsor, nstop
real(8) :: advx(0:nz+1,0:nx+1), advz(0:nz+1,0:nx+1)
real(8) :: dif1, dif2, difh, difv, diff, diffu, diffw
real(8) :: dwdx, drodx, tpp, kkk, ppp, us0, zl, dudz, drodz
real(8) :: rt, rs, ta, tb, tc, sb1, sb2, sb, rtc, rsc, tt1
real(8) :: ccs, cct, melt(0:nz+1)
real(8) :: utop, ubot
real(8) :: advy(0:nz+1,0:nx+1), um, vm, diffv, dvdz, dvdx, vtop, vbot

tax=1.3*0.0014*(wpdm)**2
pyr=2.2e-6/wpd*wpdm
intv=time/3600.
tacu=0.
do i=1, 332
  if(mod(i,4)==1)then
    tacu=tacu+12. 
    if(intv>=tacu-12. .and. intv<tacu)then
      tax=1.3*0.0014*wpd**2
	  pyr=pyr*wpd/wpdm
      exit
	endif
  elseif(mod(i,4)==2)then
    tacu=tacu+3. 
    if(intv>=tacu-3. .and. intv<tacu)then
	  wpdt=wpd-(wpd-wpdm)/3.*(intv-(tacu-3.))	
      tax=1.3*0.0014*wpdt**2
	  pyr=pyr*wpd/wpdm-(pyr*(wpd-wpdm)/wpd)/3.*(intv-(tacu-3.))
      exit
	endif	
  elseif(mod(i,4)==3)then
    tacu=tacu+17. 
    if(intv>=tacu-17. .and. intv<tacu)then
      exit
	endif
  elseif(mod(i,4)==0)then
    tacu=tacu+3. 
    if(intv>=tacu-3. .and. intv<tacu)then
	  wpdt=wpdm+(wpd-wpdm)/3.*(intv-(tacu-3.))	
      tax=1.3*0.0014*wpdt**2
	  pyr=pyr+(pyr*(wpd-wpdm)/wpd)/3.*(intv-(tacu-3.))
      exit
	endif	
  endif 
enddo
tay=0.
cc1=tay/tax
sqcd3=sqrt(0.0025)
do i=1, nx 
if(tax/=0. .or. tay/=0.) then
  w(0,i)=w(1,i)
  us(0,i)= sqrt(sqrt(tax**2+tay**2)/rho(0,i))
  ek(0,i)=us(0,i)**2/0.3 
  az(0,i)=tax*dz/sqrt(tax/rho(0,i)/sqcd3**2/sqrt(1.+cc1**2))/rho(0,i)
  taz(0,i)=az(0,i)
  if(az(0,i)< azmin) az(0,i)=azmin
  if(az(0,i)> azmax) az(0,i)=azmax 
  if(taz(0,i)< azmin) taz(0,i)=azmin
  if(taz(0,i)> azmax) taz(0,i)=azmax 
  atop=0.5*(az(1,i)+az(0,i)) 
  ep(0,i)=us(0,i)**4/atop 
  u(0,i) = u(1,i)+dz*tax/(rho(0,i)*atop)
  v(0,i) = v(1,i)+dz*tay/(rho(0,i)*atop)
  
  us0=us(0,i)
  do j=1, ni
    netwi=wi(j)+w(1,i)
    if(netwi>0.)then
      usuc=us0**2/uc2(j)
      if(usuc<1.) then
        fc(j,0,i)=fc(j,1,i)*(1.+netwi*dz/atop*usuc)
        pre(j,i)=-netwi*fc(j,1,i)*(1.-usuc)
      else
        fc(j,0,i)=fc(j,1,i)*(1.+netwi*dz/atop)
        pre(j,i)=0.
      endif
      tpre(i)=tpre(i)-pre(j,i)*dt
    else
      fc(j,0,i)=fc(j,1,i)
      pre(j,i)=0.
    endif
  enddo
else 
  u(0,i)=u(1,i)
  v(0,i)=v(1,i)  
  w(0,i)=w(1,i)  
  ek(0,i)=ek(1,i)
  ep(0,i)=ep(1,i)
  az(0,i)=az(1,i)
  taz(0,i)=taz(1,i)
  if(az(0,i)< azmin) az(0,i)=azmin
  if(az(0,i)> azmax) az(0,i)=azmax 
  if(taz(0,i)< azmin) taz(0,i)=azmin
  if(taz(0,i)> azmax) taz(0,i)=azmax   
  us(0,i)=0.
  do j=1, ni
    netwi=wi(j)+w(1,i)
    if(netwi>0.)then
      fc(j,0,i)=fc(j,1,i)
      pre(j,i)=-netwi*fc(j,1,i)
      tpre(i)=tpre(i)-pre(j,i)*dt
    else
      fc(j,0,i)=fc(j,1,i)
      pre(j,i)=0. 
    endif
  enddo
endif 

if(i==1) tpre(i-1)=tpre(i)
if(i==nx) tpre(i+1)=tpre(i)

  atop = 0.5*(taz(1,i)+taz(0,i))  
  ccs=rhoi/rho0*pyr*dz/atop
  cct=rhoi/rho0*pyr*lf/c0*dz/atop
  if((t(1,i)-(l1*s(1,i)+l2+l3*dz))<0.)then
    cct=cct*0.10
	ccs=ccs*(1.-0.10)
  endif

  dsf=1.0e-9

  if(tax/=0.)then
    s(0,i)= s(1,i)/(1.-ccs*0.50)
    if((t(1,i)-(l1*s(1,i)+l2+l3*dz))>=0.) s(0,i)=s(1,i) 
    t(0,i)= t(1,i)-cct 
  endif   
  
  tfc(0,i)=0.
  do j=1, ni
    tfc(0,i)=tfc(0,i)+fc(j,0,i)
  enddo
  rho(0,i)=rho0*(1.+bs*(s(0,i)-s0)-bt*(t(0,i)-t0)) 
  rho(0,i)=rho(0,i)+tfc(0,i)*(rhoi-rho(0,i))

  u(nz+1,i)=u(nz,i)
  v(nz+1,i)=v(nz,i)
  w(nz+1,i)=0.  
  sqcd2=sqrt(0.0025)
  tbu(i)=0.
  tbv(i)=0.
  us(nz+1,i)=0.
  az(nz+1,i)=az(nz,i)
  ek(nz+1,i)=ek(nz,i)
  ep(nz+1,i)=ep(nz,i)
  if(ek(nz+1,i)==0.) ek(nz+1,i)=ekmin
  if(ep(nz+1,i)==0.) ep(nz+1,i)=epmin
  taz(nz+1,i)=az(nz+1,i)
  if(az(nz+1,i)< azmin) az(nz+1,i)=azmin
  if(az(nz+1,i)> azmax) az(nz+1,i)=azmax 
  if(taz(nz+1,i)< azmin) taz(nz+1,i)=azmin
  if(taz(nz+1,i)> azmax) taz(nz+1,i)=azmax   
  t(nz+1,i) = t(nz,i) 
  s(nz+1,i) = s(nz,i) 
  tfc(nz+1,i)=0.
  do j=1, ni
    fc(j,nz+1,i)=fc(j,nz,i)
    tfc(nz+1,i)=tfc(nz+1,i)+fc(j,nz+1,i)
  enddo 
  rho(nz+1,i)=rho0*(1.+bs*(s(nz+1,i)-s0)-bt*(t(nz+1,i)-t0)) 
  rho(nz+1,i)=rho(nz+1,i)+tfc(nz+1,i)*(rhoi-rho(nz+1,i))
enddo

do i= 1, nz 
  if (i<=isdp) then
    w(i,0) = 0. 
    v(i,0) = 0. 
    sqcd=sqrt(0.0025)
    thv(i)=sqcd**2*v(i,1)*sqrt(1./16.*(w(i,1)+w(i,0)+w(i-1,1)+w(i-1,0))**2+v(i,1)**2)
    thw(i)=sqcd**2*w(i,1)*sqrt(1./16.*(v(i,1)+v(i,2)+v(i+1,1)+v(i+1,2))**2+w(i,1)**2)
    us(i,0)=(thw(i)**2+thv(i)**2)**0.25
    ah(i,0)=us(i,0)*sqcd*dx
    aht(i,0)=ah(i,0)
    ek(i,0)=ek(i,1)
    ep(i,0)=ep(i,1)  
    az(i,0)=az(i,1)     
    taz(i,0)=taz(i,1)
    if(az(i,0)< azmin) az(i,0)=azmin
    if(az(i,0)> azmax) az(i,0)=azmax    
    if(taz(i,0)< azmin) taz(i,0)=azmin
    if(taz(i,0)> azmax) taz(i,0)=azmax        
    
    do j=1, ni
      fc(j,i,0)=fc(j,i,1)  
    enddo    

    us0=us(i,0)  
      if((us0/sqcd)<0.05)then
        rt=1.d-3    
      else
        rt=1.d-3+gat*sqcd*(us0/sqcd-0.05)  
      endif    
      rs=0.07*rt      
      cii=ci/c0
      lff=lf/c0
      kxt=(aht(i,0)+aht(i,1))*0.5
      kxs=kxt
      sanj=rt*kxt/dx*t(i,1)/(rt+kxt/dx)
      wujx=rt**2/(rt+kxt/dx)-rt
      shiz=kxs/dx*s(i,1)/(kxs/dx+rs)
      yuel=l2+l3*i*dz-ti
      lingx=kxs/dx+rs
      ta=l1*cii*rs**2/lingx-l1*rs*cii-l1*wujx
      tb=rs**2/lingx*cii*yuel+rs*shiz*l1*cii-rs*lff-rs*cii*yuel+rs**2/lingx*lff-sanj-(l2+l3*i*dz)*wujx
      tc=rs*shiz*lff+rs*shiz*cii*yuel
      sb1=(-tb+sqrt(tb**2-4.*ta*tc))/2./ta
      sb2=(-tb-sqrt(tb**2-4.*ta*tc))/2./ta
      if(sb1>0. .and. sb1<50.) sb=sb1
      if(sb2>0. .and. sb2<50.) sb=sb2
      tb=l1*sb+l2+l3*i*dz
      u(i,0)=rs/sb*(rs*sb+kxs*s(i,1)/dx)/(kxs/dx+rs)-rs
      melt(i)=u(i,0)         
      t(i,0)=kxt/dx*t(i,1)/(rt+kxt/dx)+rt/(rt+kxt/dx)*tb
      s(i,0)=(rs*sb+kxs/dx*s(i,1))/(kxs/dx+rs)
      if(ad==1.) tfre(i)=tfre(i)+u(i,0)*dt*rho0/rhoi 
  elseif (i>isdp .and. i<=isdp2) then
    if(n==1)then
      tt0=t(i,0)
      ss0=s(i,0)
    endif
    if(i<=isdp+isd)then
      if(i==isdp+1)then 
        u(i,0)=ad*uini  
        t(i,0)=tt0+ad*(tmix-tt0) 
        s(i,0)=ss0+ad*(smix-ss0) 
        do j=1, ni
          fc(j,i,0)=fc(j,i,1) 
        enddo
      else
        u(i,0)=u(i-1,0)          
        t(i,0)=t(i-1,0)
        s(i,0)=s(i-1,0)
        do j=1, ni
          fc(j,i,0)=fc(j,i-1,0)
        enddo
      endif
    else
      u(i,0)=u(i,1) 
      t(i,0)=t(i,1)
      s(i,0)=s(i,1)
      do j=1, ni
        fc(j,i,0)=fc(j,i,1)  
      enddo         
	endif 
    v(i,0)=v(i,1)    
    w(i,0)=w(i,1)	
    ek(i,0)=ek(i,1)
    ep(i,0)=ep(i,1)
    az(i,0)=az(i,1)
    taz(i,0)=taz(i,1)
    if(az(i,0)< azmin) az(i,0)=azmin
    if(az(i,0)> azmax) az(i,0)=azmax    
    if(taz(i,0)< azmin) taz(i,0)=azmin
    if(taz(i,0)> azmax) taz(i,0)=azmax 
  else
    u(i,0)=u(i,1)
    v(i,0)=v(i,1)    
    w(i,0)=w(i,1)    
    t(i,0)=t(i,1)
    s(i,0)=s(i,1)
    do j=1, ni
      fc(j,i,0)=fc(j,i,1)  
    enddo    
    ek(i,0)=ek(i,1)
    ep(i,0)=ep(i,1)
    az(i,0)=az(i,1)
    taz(i,0)=taz(i,1)
    if(az(i,0)< azmin) az(i,0)=azmin
    if(az(i,0)> azmax) az(i,0)=azmax    
    if(taz(i,0)< azmin) taz(i,0)=azmin
    if(taz(i,0)> azmax) taz(i,0)=azmax 
  endif
  tfc(i,0)=0.
  do j=1, ni
    tfc(i,0)=tfc(i,0)+fc(j,i,0)
  enddo 
  rho(i,0)=rho0*(1.+bs*(s(i,0)-s0)-bt*(t(i,0)-t0)) 
  rho(i,0)=rho(i,0)+tfc(i,0)*(rhoi-rho(i,0))  

    u(i,nx+1)=u(i,nx) 
    v(i,nx+1)=v(i,nx)    
    w(i,nx+1)=w(i,nx)
    t(i,nx+1)=t(i,nx)
    s(i,nx+1)=s(i,nx)
    do j=1, ni
      fc(j,i,nx+1)=fc(j,i,nx)  
    enddo    
    ek(i,nx+1)=ek(i,nx)
    ep(i,nx+1)=ep(i,nx)
    az(i,nx+1)=az(i,nx)
    taz(i,nx+1)=taz(i,nx)
    if(az(i,nx+1)< azmin) az(i,nx+1)=azmin
    if(az(i,nx+1)> azmax) az(i,nx+1)=azmax  
    if(taz(i,nx+1)< azmin) taz(i,nx+1)=azmin
    if(taz(i,nx+1)> azmax) taz(i,nx+1)=azmax 
    tfc(i,nx+1)=0.
    do j=1, ni
      fc(j,i,nx+1)=fc(j,i,nx)
      tfc(i,nx+1)=tfc(i,nx+1)+fc(j,i,nx+1)
    enddo 
    rho(i,nx+1)=rho0*(1.+bs*(s(i,nx+1)-s0)-bt*(t(i,nx+1)-t0)) 
    rho(i,nx+1)=rho(i,nx+1)+tfc(i,nx+1)*(rhoi-rho(i,nx+1))      
enddo

IF(MOD(n,nout)==0)THEN
  
  ntem=n/nout 
  write( cTemp,'(i3)' ) ntem

  open ( 1, file = 'uw' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,u,w"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15,x,f20.15)') (k-1)*dx, -(i-1)*dz, u(i,k), w(i,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 't' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,t"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz 
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz, t(i,k)
    enddo
  END DO  
  close(1)     
  
  open ( 1, file = 's' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,s"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz 
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz, s(i,k)
    enddo
  END DO  
  close(1)   
  
  open ( 1, file = 'tsc' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,tsc"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz,  t(i,k)-(l1*s(i,k)+l2+l3*i*dz)
    enddo
  END DO  
  close(1)  
  
  open ( 1, file = 'rho' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,rho"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz 
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz, rho(i,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'tfc' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,tfc"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz 
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz, tfc(i,k)
    enddo
  END DO  
  close(1)     

  open ( 1, file = 'sfc' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,sfc"
  write(1,*) "zone i=", nx, "j=", nz, "f=point"
  DO i = 1, nz 
    do k= 1, nx
      WRITE(1,'(f7.1,x,f7.1,x,f20.15)') (k-1)*dx, -(i-1)*dz, sfc(i,k)
    enddo
  END DO  
  close(1)      
 
ENDIF

DO i = 0,nz+1
DO k = 0,nx+1
  CuP(i,k) = 0.5*(u(i,k)+abs(u(i,k)))*dt/dx
  CuN(i,k) = 0.5*(u(i,k)-abs(u(i,k)))*dt/dx
  CwP(i,k) = 0.5*(w(i,k)+abs(w(i,k)))*dt/dz
  CwN(i,k) = 0.5*(w(i,k)-abs(w(i,k)))*dt/dz
END DO
END DO

DO i = 1,nz
DO k = 1,nx
 divv(i,k) = (u(i,k)-u(i,k-1))/dx+(w(i,k)-w(i+1,k))/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = ek(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div = dt*B(i,k)*divv(i,k)
 azt(i,k)=0.5*(az(i,k)+az(i-1,k))
 IF(WET(i-1,k)) dif1 = azt(i,k)*(ek(i-1,k)-ek(i,k))/dz
 azb(i,k)=0.5*(az(i,k)+az(i+1,k))
 dif2 = azb(i,k)*(ek(i,k)-ek(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diff = dt*difz  
 dudz = (u(i-1,k)-u(i+1,k))/(2.*dz) 
 dvdz = (v(i-1,k)-v(i+1,k))/(2.*dz)
 shear(i,k)=dudz**2+ dvdz**2
 ps(i,k) = shear(i,k)*az(i,k)
 nn22(i,k)=gv/rho0*(rho(i-1,k)-rho(i+1,k))/(2.*dz)
 if(az(i,k)==azmin)then
   taz(i,k)=azmin  
 else    
   if(shear(i,k)==0.) then  
     taz(i,k)=azmin  
   else  
     rich(i,k)=-nn22(i,k)/shear(i,k)
     if(rich(i,k)<=0.)then
       taz(i,k)= az(i,k)/0.7 
     else
     taz(i,k)=az(i,k)/(0.7*(1.+10.*rich(i,k))**(-0.5)/(1.+10./3.*rich(i,k))**(-1.5))
     endif
     if(taz(i,k)<azmin) taz(i,k)=azmin   
   endif 
 endif
 pb(i,k)=gv*(-bt*taz(i,k)*(t(i-1,k)-t(i+1,k))/(2.*dz)+bs*taz(i,k)*(s(i-1,k)-s(i+1,k))/(2.*dz))
 pb(i,k)=pb(i,k)+az(i,k)*gv*(rhoi-rho0)/rho0*(tfc(i-1,k)-tfc(i+1,k))/(2.*dz)
 kkk=ek(i,k)
 ekn(i,k)=ek(i,k)+BN(i,k)+div+diff/1.4+dt*(ps(i,k)+pb(i,k)-ep(i,k))
 if(ekn(i,k)<0.) ekn(i,k)=kkk

END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = ep(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div = dt*B(i,k)*divv(i,k)
 IF(WET(i-1,k)) dif1 = azt(i,k)*(ep(i-1,k)-ep(i,k))/dz
 dif2 = azb(i,k)*(ep(i,k)-ep(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diff = dt*difz  
 ppp=ep(i,k)
 epn(i,k)=ep(i,k)+BN(i,k)+div+diff/1.3+dt*ep(i,k)/ek(i,k)*(1.44*ps(i,k)+0.8*pb(i,k)-1.92*ep(i,k))
 if(epn(i,k)<0.) epn(i,k)=ppp

END DO
END DO

DO i = 1,nz
DO k = 1,nx
 ek(i,k)=ekn(i,k)   
 ep(i,k)=epn(i,k)
END DO
END DO

DO i = 1, nz
DO k = 1, nx
  az(i,k)=0.09*ek(i,k)**2/ep(i,k)
  if(az(i,k)< azmin) az(i,k)=azmin
  if(az(i,k)> azmax) az(i,k)=azmax  
  if(ek(i,k)<=ekmin) az(i,k)=azmin
END DO
END DO 

ttfc=0.
sfc=0.
if (cmin==0.) goto 888 

do i=1,nz
do k=1,nx
    
  ttfc(i,k)=0.
  sfc(i,k)=0.
  senu(1,i,k)=0.    

  do j=1,ni  

    rtc=kt*nu(j)/ri(j)*1000.  
    rsc=ks*sh(j)/ri(j)*1000.  
    ta=l1*rtc*c0
    tb=(l2+l3*i*dz-t(i,k))*c0*rtc-rsc*lf    
    tc=rsc*lf*s(i,k)
    sb1=(-tb+sqrt(tb**2-4.*ta*tc))/2./ta
    sb2=(-tb-sqrt(tb**2-4.*ta*tc))/2./ta
    if(sb1>0. .and. sb1<50.) sb=sb1
    if(sb2>0. .and. sb2<50.) sb=sb2    
    wf0(j,i,k)=-(1.-fc(j,i,k))*rsc*(s(i,k)-sb)*2.*fc(j,i,k)/ri(j)/sb   
    if(wf0(j,i,k)<0.) then
      wf0(j,i,k)=wf0(j,i,k)*(1.+0.5/ar)  
      if(abs(wf0(j,i,k))*dt>fc(j,i,k))  wf0(j,i,k)=-fc(j,i,k)/dt
    endif   
    ttfc(i,k)=ttfc(i,k)-(l1*sb+l2+l3*i*dz-t(i,k)-lf/c0)*wf0(j,i,k)
    sfc(i,k)=sfc(i,k)+s(i,k)*wf0(j,i,k) 
    wt(j)=sqrt(wi(j)**2+4.*ep(i,k)*re(j)**2/15./1.95e-6)
    if(j>1) then
      senu(j,i,k)=-pi*nni*wt(j)/re(j)*re(1)**3*fc(j,i,k)
      senu(1,i,k)=senu(1,i,k)-senu(j,i,k)  
    endif
  enddo
  
  do j=1,ni
    if (wf0(j,i,k)>=0.) then  
      if(j==1) then
        wf(j,i,k)=wf0(j,i,k)/dv(j)*vi(j)
      elseif(j>1 .and. j<ni) then
        if(wf0(j-1,i,k)<0.) then
          wf(j,i,k)=wf0(j,i,k)/dv(j)*vi(j)
        else
          wf(j,i,k)=-(wf0(j-1,i,k)/dv(j-1)-wf0(j,i,k)/dv(j))*vi(j)
        endif
      elseif(j==ni) then
        if(wf0(j-1,i,k)<0.) then
          wf(j,i,k)=-wf0(j,i,k) 
        else
          wf(j,i,k)=-wf0(j,i,k)-wf0(j-1,i,k)/dv(j-1)*vi(j)
        endif
      endif
    else
      if(j==1) then
        if(wf0(j+1,i,k)>0.) then
          wf(j,i,k)=-wf0(j,i,k)  
        else
          wf(j,i,k)=wf0(j+1,i,k)/dv(j)*vi(j)-wf0(j,i,k)
        endif
      elseif(j>1 .and. j<ni) then
        if(wf0(j+1,i,k)>0.) then
          wf(j,i,k)=-wf0(j,i,k)/dv(j-1)*vi(j)
        else 
          wf(j,i,k)=(wf0(j+1,i,k)/dv(j)-wf0(j,i,k)/dv(j-1))*vi(j)
        endif  
      elseif(j==ni) then
        wf(j,i,k)=-wf0(j,i,k)/dv(j-1)*vi(j) 
      endif      
    endif    
  enddo
enddo
enddo

do j=1, ni
  DO i = 0,nz+1
  DO k = 0,nx+1
    B(i,k) = fc(j,i,k)
  END DO
  END DO

  CALL advect

  DO i = 1,nz
  DO k = 1,nx
    div = dt*B(i,k)*divv(i,k)
    IF(WET(i,k+1)) dif1 = 0.5*(ah(i,k)+ah(i,k+1))*(fc(j,i,k+1)-fc(j,i,k))/dx
    dif2=0.
    IF(WET(i,k-1)) dif2 = 0.5*(ah(i,k)+ah(i,k-1))*(fc(j,i,k)-fc(j,i,k-1))/dx
    difh = (dif1-dif2)/dx
    IF(WET(i-1,k)) dif1 = azt(i,k)*(fc(j,i-1,k)-fc(j,i,k))/dz
    dif2=0.
    IF(WET(i+1,k)) dif2 = azb(i,k)*(fc(j,i,k)-fc(j,i+1,k))/dz
    difz = (dif1-dif2)/dz
    diff = dt*(difh+difz)  
    fcn(j,i,k)= fc(j,i,k)+BN(i,k)+div+diff -dt*wi(j)*(fc(j,i-1,k)-fc(j,i+1,k))/2./dz -rho0/rhoi*wf(j,i,k)*dt +senu(j,i,k)*dt
    if(j==1 .and. i==1 .and. (t(1,k)-(l1*s(1,k)+l2+l3*dz))<0.)  fcn(j,i,k)= fcn(j,i,k) +dsf*dt  
    if (fcn(j,i,k)<0.) fcn(j,i,k)= 0.
    
  END DO
  END DO

enddo

DO i = 1,nz
DO k = 1,nx
  tfc(i,k)=0.	
  do j=1,ni	
    tfc(i,k)=tfc(i,k)+fcn(j,i,k)
  enddo
  if(tfc(i,k)>1.e-3)then
    do j=1,ni	
      fcn(j,i,k)=fc(j,i,k)
    enddo  
  endif
  do j=1,ni	
    fc(j,i,k)=fcn(j,i,k)
  enddo   
END DO
END DO

DO i = 0, nz+1
DO k = 0, nx+1
  tfc(i,k)=0.
  do j=1,ni
    tfc(i,k)=tfc(i,k)+fc(j,i,k)  
  enddo
END DO
END DO

888 DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = t(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div = dt*B(i,k)*divv(i,k)
 IF(WET(i,k+1)) dif1 = 0.5*(aht(i,k)+aht(i,k+1))*(t(i,k+1)-t(i,k))/dx
 dif2 = 0.5*(aht(i,k)+aht(i,k-1))*(t(i,k)-t(i,k-1))/dx
 difh = (dif1-dif2)/dx
 aztt(i,k) = 0.5*(taz(i,k)+taz(i-1,k))
 IF(WET(i-1,k)) dif1 = aztt(i,k)*(t(i-1,k)-t(i,k))/dz
 dif2=0.
 azbb(i,k) = 0.5*(taz(i,k)+taz(i+1,k))
 IF(WET(i+1,k)) dif2 = azbb(i,k)*(t(i,k)-t(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 tn(i,k)= t(i,k)+BN(i,k)+div+diff +ttfc(i,k)*dt
END DO
END DO

DO i = 1,nz
DO k = 1,nx
 t(i,k)=tn(i,k)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = s(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
 div = dt*B(i,k)*divv(i,k)
 IF(WET(i,k+1)) dif1 = 0.5*(aht(i,k)+aht(i,k+1))*(s(i,k+1)-s(i,k))/dx
 dif2 = 0.5*(aht(i,k)+aht(i,k-1))*(s(i,k)-s(i,k-1))/dx
 difh = (dif1-dif2)/dx
 IF(WET(i-1,k)) dif1 = aztt(i,k)*(s(i-1,k)-s(i,k))/dz
 dif2=0.
 IF(WET(i+1,k)) dif2 = azbb(i,k)*(s(i,k)-s(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diff = dt*(difh+difz)  
 ssn(i,k)= s(i,k)+BN(i,k)+div+diff +sfc(i,k)*dt
END DO
END DO

DO i = 1,nz
DO k = 1,nx
 s(i,k)=ssn(i,k)
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  rho(i,k) = rho0*(1.+bs*(s(i,k)-s0)-bt*(t(i,k)-t0)) 
  rho(i,k) = rho(i,k)+tfc(i,k)*(rhoi-rho(i,k))
END DO
END DO

DO k = 0,nx+1
  p(0,k) = 0.0  
  DO i = 1,nz+1
    p(i,k) = p(i-1,k) + (0.5*(rho(i-1,k)+rho(i,k)))*gv*dz
  END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = v(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
  div1 = (u(i,k)-u(i,k-1))/dx
  div2 = (w(i,k)-w(i+1,k))/dz
  div = dt*B(i,k)*(div1+div2)  
  advy(i,k)= BN(i,k)+div
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx
  CuP(i,k) = 0.25*(u(i,k)+u(i,k+1)+abs(u(i,k))+abs(u(i,k+1)))*dt/dx
  CuN(i,k) = 0.25*(u(i,k)+u(i,k+1)-abs(u(i,k))-abs(u(i,k+1)))*dt/dx
  CwP(i,k) = 0.25*(w(i,k)+w(i,k+1)+abs(w(i,k))+abs(w(i,k+1)))*dt/dz
  CwN(i,k) = 0.25*(w(i,k)+w(i,k+1)-abs(w(i,k))-abs(w(i,k+1)))*dt/dz
END DO
END DO
  
DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = u(i,k)
END DO
END DO
  
CALL advect  

DO i = 1,nz
DO k = 1,nx
  div1 = 0.5*(u(i,k+1)-u(i,k-1))/dx
  div2 = 0.5*(w(i,k)+w(i,k+1)-w(i+1,k)-w(i+1,k+1))/dz
  div = dt*B(i,k)*(div1+div2)  
  advx(i,k)= BN(i,k)+div
END DO
END DO

DO i = 0,nz
DO k = 0,nx+1
  CuP(i,k) = 0.25*(u(i,k)+u(i-1,k)+abs(u(i,k))+abs(u(i-1,k)))*dt/dx
  CuN(i,k) = 0.25*(u(i,k)+u(i-1,k)-abs(u(i,k))-abs(u(i-1,k)))*dt/dx
  CwP(i,k) = 0.25*(w(i,k)+w(i-1,k)+abs(w(i,k))+abs(w(i-1,k)))*dt/dz
  CwN(i,k) = 0.25*(w(i,k)+w(i-1,k)-abs(w(i,k))-abs(w(i-1,k)))*dt/dz
END DO
END DO

DO i = 0,nz+1
DO k = 0,nx+1
  B(i,k) = w(i,k)
END DO
END DO

CALL advect

DO i = 1,nz
DO k = 1,nx
  div1 = 0.5*(u(i,k)+u(i-1,k)-u(i,k-1)-u(i-1,k-1))/dx
  div2 = 0.5*(w(i-1,k)-w(i+1,k))/dz
  div = dt*B(i,k)*(div1+div2)  
  advz(i,k)= BN(i,k) + div
END DO
END DO

DO i = 1,nz
DO k = 1,nx

 IF(WET(i,k+1)) dif1 = ah(i,k+1)*(u(i,k+1)-u(i,k))/dx
 dif2=0.
 IF(WET(i,k-1)) dif2 = ah(i,k)*(u(i,k)-u(i,k-1))/dx
 difh = (dif1-dif2)/dx 
 atop=0.25*(az(i,k)+az(i-1,k)+az(i,k+1)+az(i-1,k+1))
 IF(WET(i-1,k)) dif1 = atop*(u(i-1,k)-u(i,k))/dz
 if(i==1) dif1=tax/rho(0,k)
 abot=0.25*(az(i,k)+az(i+1,k)+az(i,k+1)+az(i+1,k+1))
 dif2 = abot*(u(i,k)-u(i+1,k))/dz
 if(i==nz) dif2=tbu(k)
 difz = (dif1-dif2)/dz
 diffu = dt*(difh+difz)  

 IF(WET(i,k+1)) dif1 = 0.5*(ah(i,k)+ah(i,k+1))*(v(i,k+1)-v(i,k))/dx
 IF(WET(i,k-1)) dif2 = 0.5*(ah(i,k)+ah(i,k-1))*(v(i,k)-v(i,k-1))/dx
 if(k==1 .and. i<=isdp) dif2 = thv(i)
 difh = (dif1-dif2)/dx
 atop=0.5*(az(i,k)+az(i-1,k))
 IF(WET(i-1,k)) dif1 = atop*(v(i-1,k)-v(i,k))/dz
 if(i==1) dif1=tay/rho(0,k)
 abot=0.5*(az(i,k)+az(i+1,k))
 dif2 = abot*(v(i,k)-v(i+1,k))/dz
 if(i==nz) dif2=tbv(k) 
 difz = (dif1-dif2)/dz
 diffv = dt*(difh+difz)

 IF(WET(i,k+1)) dif1 = 0.5*(ah(i,k)+ah(i,k+1))*(w(i,k+1)-w(i,k))/dx
 IF(WET(i,k-1)) dif2 = 0.5*(ah(i,k)+ah(i,k-1))*(w(i,k)-w(i,k-1))/dx
 if(k==1 .and. i<=isdp) dif2 = thw(i)
 difh = (dif1-dif2)/dx
 IF(WET(i-1,k)) dif1 = az(i-1,k)*(w(i-1,k)-w(i,k))/dz
 dif2=0.
 IF(WET(i+1,k)) dif2 = az(i,k)*(w(i,k)-w(i+1,k))/dz
 difz = (dif1-dif2)/dz
 diffw = dt*(difh+difz)  
 
 IF(wet(i,k))THEN
     
   um = 0.5*(u(i,k)+u(i,k-1))
   vm = 0.5*(v(i,k)+v(i,k+1))     
     
   pressx = -drdx*(q(i,k+1)-q(i,k))-drdx*(p(i,k+1)-p(i,k))
   IF(wet(i,k+1)) ustar(i,k) = cos(alpha)*u(i,k) + sin(alpha)*vm + dt*pressx + advx(i,k) + diffu
   
   vstar(i,k) = cos(alpha)*v(i,k) - sin(alpha)*um + advy(i,k) + diffv
   
   pressz = -drdz*(q(i-1,k)-q(i,k))
   IF(wet(i-1,k)) wstar(i,k) = w(i,k) + dt*pressz + advz(i,k) + diffw
 ENDIF
END DO
END DO

DO i = 1,nz
 if(i<=isdp) then
   ustar(i,0)=melt(i)
   vstar(i,0)=0.
   wstar(i,0)=0.
 elseif(i>isdp .and. i<=isdp2) then
    if(i<=isdp+isd)then
      if(i==isdp+1)then 
        ustar(i,0)=ad*uini        
      else
        ustar(i,0)=ustar(i-1,0)         
      endif
      vstar(i,0) = vstar(i,1)
      wstar(i,0) = wstar(i,1) 
    else
      ustar(i,0)=ustar(i,1)
      vstar(i,0)=vstar(i,1)    
      wstar(i,0)=wstar(i,1)          
    endif    
 else
   ustar(i,0) = ustar(i,1) 
   vstar(i,0) = vstar(i,1)
   wstar(i,0) = wstar(i,1)
 endif

 ustar(i,nx+1) = ustar(i,nx)
 vstar(i,nx+1) = vstar(i,nx) 
 wstar(i,nx+1) = wstar(i,nx) 
END DO

DO i = 1,nz
DO k = 1,nx
 qstar(i,k) = -1.*rho0/dt*(  &
 &  (ustar(i,k)-ustar(i,k-1))*dz + &
 &  (wstar(i,k)-wstar(i+1,k))*dx )
END DO
END DO

nstop = 8000

!*****************
DO nsor = 1,nstop
!*****************

perr = 0.0

DO i = 1,nz
DO k = 1,nx

 IF(wet(i,k))THEN     
 q1 = dq(i,k)
 term1 = qstar(i,k) + & 
  &      at(i,k)*dq(i-1,k) + ab(i,k)*dq(i+1,k) + & 
  &      aw(i,k)*dq(i,k-1) + ae(i,k)*dq(i,k+1)
 q2 = (1.0-omega)*q1 + omega*term1/atot(i,k) 
 dq(i,k) = q2
 perr = MAX(ABS(q2-q1),perr)
 ENDIF
 
END DO
END DO

DO i = 1,nz
DO k = 1,nx
  IF(wet(i,k))THEN
    pressx = -drdx*(dq(i,k+1)-dq(i,k))
   IF(wet(i,k+1)) un(i,k) = ustar(i,k) + dt*pressx
   vn(i,k) = vstar(i,k)
   pressz = -drdz*(dq(i-1,k)-dq(i,k))
   IF(wet(i-1,k)) wn(i,k) = wstar(i,k) + dt*pressz
  END IF
END DO
END DO

DO k = 1,nx
  usum(k) = 0.
  DO i = 1,nz 
    usum(k) = usum(k) + dz*un(i,k)
  END DO
END DO

usum(0) = 0.
usum(nx+1) = usum(nx)
DO i = 1,nz
  if(i<=isdp) then
    usum(0)=usum(0)+dz*melt(i)
  elseif(i>isdp .and. i<=isdp2) then
    if(i<=isdp+isd)then
      usum(0) = usum(0) +ad*uini*dz  
    else
      usum(0) = usum(0) + dz*un(i,1) 
    endif     
  else
    usum(0) = usum(0) + dz*un(i,1)
  endif
END DO

DO k = 1,nx
 dq(0,k) = -dt*rho0*gv*(usum(k)-usum(k-1))/dx
END DO

IF(perr <= peps)THEN
  nstop = nsor
  GOTO 33
END IF

!********************
END DO
!********************

GOTO 34
33 WRITE(*,*) "No. of Interactions =>", nstop

34 CONTINUE

DO i = 1,nz
DO k = 1,nx
  q(i,k) = q(i,k)+dq(i,k)
  u(i,k) = un(i,k)
  v(i,k) = vn(i,k)
  w(i,k) = wn(i,k)
END DO
END DO

DO k = 1,nx
  q(0,k) = q(0,k)+dq(0,k)
END DO

DO i = 0,nz
  q(i,0)=q(i,1)
  q(i,nx+1) = q(i,nx)
END DO


RETURN
END SUBROUTINE dyn

SUBROUTINE advect
implicit none
real(8) :: RxP(0:nz+1,0:nx+1), RxN(0:nz+1,0:nx+1)
real(8) :: RzP(0:nz+1,0:nx+1), RzN(0:nz+1,0:nx+1)
real(8) :: dB, term1, term2, term3, term4
real(8) :: BwP, BwN, BeP, BeN, BbP, BbN, BtP, BtN 

DO i = 0,nz+1
DO k = 0,nx+1
  RxP(i,k) = 0.0
  RxN(i,k) = 0.0
  RzP(i,k) = 0.0
  RzN(i,k) = 0.0
END DO
END DO

DO i = 1,nz
DO k = 1,nx
  dB =  B(i,k+1)-B(i,k)
  IF(ABS(dB) > 0.0) RxP(i,k) = (B(i,k)-B(i,k-1))/dB
  dB =  B(i-1,k)-B(i,k)
  IF(ABS(dB) > 0.0) RzP(i,k) = (B(i,k)-B(i+1,k))/dB
END DO
END DO

DO i = 1,nz
DO k = 0,nx-1
  dB =  B(i,k+1)-B(i,k)
  IF(ABS(dB) > 0.0) RxN(i,k) = (B(i,k+2)-B(i,k+1))/dB
END DO
END DO

DO i = 2,nz+1
DO k = 1,nx
  dB =  B(i-1,k)-B(i,k)
  IF(ABS(dB) > 0.0) RzN(i,k) = (B(i-2,k)-B(i-1,k))/dB
END DO
END DO   

DO i = 1,nz
DO k = 1,nx

term1 = (1.0-CuP(i,k-1))*(B(i,k)-B(i,k-1))
BwP = B(i,k-1)+0.5*PSI(RxP(i,k-1))*term1

term1 = (1.0+CuN(i,k-1))*(B(i,k)-B(i,k-1))
BwN = B(i,k)-0.5*PSI(RxN(i,k-1))*term1

term1 = (1.0-CuP(i,k))*(B(i,k+1)-B(i,k))
BeP = B(i,k)+0.5*PSI(RxP(i,k))*term1

term1 = (1.0+CuN(i,k))*(B(i,k+1)-B(i,k))  
BeN = B(i,k+1)-0.5*PSI(RxN(i,k))*term1

term1 = (1.0-CwP(i+1,k))*(B(i,k)-B(i+1,k))
BbP = B(i+1,k)+0.5*PSI(RzP(i+1,k))*term1

term1 = (1.0+CwN(i+1,k))*(B(i,k)-B(i+1,k))  
BbN = B(i,k)-0.5*PSI(RzN(i+1,k))*term1

term1 = (1.0-CwP(i,k))*(B(i-1,k)-B(i,k)) 
BtP = B(i,k)+0.5*PSI(RzP(i,k))*term1

term1 = (1.0+CwN(i,k))*(B(i-1,k)-B(i,k)) 
BtN = B(i-1,k)-0.5*PSI(RzN(i,k))*term1

term1 = CuP(i,k-1)*BwP+CuN(i,k-1)*BwN
term2 = CuP(i,k)*BeP+CuN(i,k)*BeN
term3 = CwP(i+1,k)*BbP+CwN(i+1,k)*BbN
term4 = CwP(i,k)*BtP+CwN(i,k)*BtN

BN(i,k) = term1-term2+term3-term4

END DO
END DO

RETURN

END SUBROUTINE advect

REAL FUNCTION psi(r)
implicit none
real(8), INTENT(IN) :: r  
real(8) :: term1, term2, term3

term1 = MIN(2.0*r,1.0)
term2 = MIN(r,2.0)
term3 = MAX(term1,term2)
psi = MAX(term3,0.0)

RETURN

END FUNCTION psi 

END MODULE sub
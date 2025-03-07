PROGRAM slice
USE param
USE sub

!**********
CALL INIT  
!**********

ntot = int(121.*24.*3600./dt)+1 
time = 0.0
ntt=4.
nout = int(ntt*3600./dt)+1 

DO n = 1, ntot
  time = time + dt
  write(6,*)"time (hours)", time/(3600.)
  ad = 1.
  isd=min((isdp2-isdp), int((isdp2-isdp)*ad))
  isd=max(isd, 1)
  CALL dyn

END DO

END PROGRAM slice

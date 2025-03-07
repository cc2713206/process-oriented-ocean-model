MODULE param
implicit none
INTEGER(4), PARAMETER :: nx = 201
INTEGER(4), PARAMETER :: nz = 81
real, parameter :: gv = 9.81 
real(8) :: tax
real(8) :: tay
integer, parameter :: ni=9
real, parameter :: dx = 50.
real, parameter :: dz = 10.
real, parameter :: dt = 1.6
real(8), PARAMETER :: tta = -75. 
real(8), parameter :: pi=3.14159265359
real(8), parameter :: om=7.292e-5
integer, parameter :: isdp=29
integer, parameter :: isdp2=33
real, parameter :: ar=1./50.
real, parameter :: rho0=1030.
real, parameter :: rhoi=920.
real, parameter :: t0=-2.
real, parameter :: s0=34.5
real, parameter :: bt=3.87e-5
real, parameter :: bs=7.86e-4
real, parameter :: cmin=1.e-7
real(8) :: dsf
real, parameter :: uini = 0.02 
real, parameter :: fini = 2.e-5 
real(8) :: tmix
real, parameter :: smix=34.7
real, parameter :: gat=2.2e-2
real, parameter :: gas=6.2e-4
real, parameter :: l1=-0.0573
real, parameter :: l2=0.0832
real, parameter :: l3=-7.61e-4
real, parameter :: c0=3974.
real, parameter :: lf=3.35e5
real, parameter :: r1=0.01
real, parameter :: r2=0.02
real, parameter :: r3=0.05
real, parameter :: r4=0.1
real, parameter :: r5=0.2
real, parameter :: r6=0.3
real, parameter :: r7=0.5
real, parameter :: r8=0.8
real, parameter :: r9=1.0
real, parameter :: kt=1.4e-7
real, parameter :: ks=8.d-10
real, parameter :: theta=0.05 
real, parameter :: tsea=0.015 
real, parameter :: azmin=1.e-5
real, parameter :: ekmin=1.e-9 
real, parameter :: epmin=1.e-12 
real, parameter :: azmax=0.5
real, parameter :: nni=1000. 
real(8) :: ustar(0:nz+1,0:nx+1), wstar(0:nz+1,0:nx+1), qstar(nz,nx)
real, parameter :: omega=1.5 
real, parameter :: peps=1.e-3 
real(8) :: at(nz,nx), ab(nz,nx) 
real(8) :: ae(nz,nx), aw(nz,nx)
real(8) :: atot(nz,nx) 
real(8) :: drdz, drdx 
INTEGER :: n,i,j,k
real(8) :: tfe(5000, 0:nz+1), tpe(5000, 0:nx+1)
real(8) :: time
real(8) :: rho(0:nz+1,0:nx+1),rhon(0:nz+1,0:nx+1), difz
real(8) :: p(0:nz+1,0:nx+1) 
real(8) :: dq(0:nz+1,0:nx+1) 
real(8) :: q(0:nz+1,0:nx+1) 
real(8) :: u(0:nz+1,0:nx+1), un(0:nz+1,0:nx+1) 
real(8) :: w(0:nz+1,0:nx+1), wn(0:nz+1,0:nx+1) 
real(8) :: usum(0:nx+1) 
LOGICAL :: wet(0:nz+1,0:nx+1) 
real(8) :: CuP(0:nz+1,0:nx+1), CuN(0:nz+1,0:nx+1)
real(8) :: CwP(0:nz+1,0:nx+1), CwN(0:nz+1,0:nx+1)
real(8) :: B(0:nz+1,0:nx+1), BN(0:nz+1,0:nx+1)
real(8) :: wi(ni), nu(ni), sh(ni), ri(ni), re(ni), uc2(ni), vi(ni), wt(ni), dv(ni)
real(8) :: tfc(0:nz+1,0:nx+1), fc(ni,0:nz+1,0:nx+1), t(0:nz+1,0:nx+1), s(0:nz+1,0:nx+1), tsc(0:nz+1,0:nx+1)
real(8) :: tfre(0:nz+1), tpre(0:nx+1), hsc(0:nx+1), pre(ni,0:nx+1), sea(ni,0:nx+1)
real(8) :: az(0:nz+1,0:nx+1), us(0:nz+1,0:nx+1), ek(0:nz+1,0:nx+1), ep(0:nz+1,0:nx+1)
real(8) :: div, div1, div2, divv(0:nz+1,0:nx+1)
real(8) :: atop, abot, ps(0:nz+1,0:nx+1), pb(0:nz+1,0:nx+1)
real(8) :: ttfc(0:nz+1,0:nx+1), sfc(0:nz+1,0:nx+1), senu(ni,0:nz+1,0:nx+1), wf0(ni,0:nz+1,0:nx+1), wf(ni,0:nz+1,0:nx+1)
real(8) :: dzz
INTEGER :: ntot, nout
real :: ntt
integer :: ntem
character( len = 3 ) :: cTemp
real(8) :: ekn(0:nz+1,0:nx+1), epn(0:nz+1,0:nx+1), fcn(ni,0:nz+1,0:nx+1), tn(0:nz+1,0:nx+1),ssn(0:nz+1,0:nx+1)
real(8) :: v(0:nz+1,0:nx+1), vn(0:nz+1,0:nx+1), vstar(0:nz+1,0:nx+1), azt(0:nz+1,0:nx+1), azb(0:nz+1,0:nx+1)
real(8) :: alpha
real(8) :: sqcd
integer :: isd
real, parameter :: ci=2009.
real, parameter :: ti=-25.
real (8) :: lzi,  cii, lff, kxt, kxs,sanj, wujx, shiz,yuel,lingx, usuc, pyr, sqcd2, cc1, sqcd3, ad, eddyx, eddyy, stab, wurz, c2, netwi, tt0, ss0 
real (8) :: tbu(0:nx+1), tbv(0:nx+1), thv(0:nz+1), thw(0:nz+1), taz(0:nz+1,0:nx+1), ah(0:nz+1,0:nx+1), aht(0:nz+1,0:nx+1), shear(0:nz+1,0:nx+1), rich(0:nz+1,0:nx+1), nn22(0:nz+1,0:nx+1), aztt(0:nz+1,0:nx+1), azbb(0:nz+1,0:nx+1)
LOGICAL :: dry(0:nz+1,0:nx+1)
real (8) :: intv, tacu, wpdt
real, parameter :: wpd=33.
real, parameter :: wpdm=4.
END MODULE param
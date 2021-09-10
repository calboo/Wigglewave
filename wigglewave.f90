!                                          
!   \  X  / o  _   _   )  _         _      _  
!    \/ \/  ( (_( (_( (  )_) )_)_) (_( \) )_) 
!               _)  _)  (_               (_   
!

! Copyright (C) 2020      Callum Boocok <c.boocock@qmul.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details:
! <http://www.gnu.org/licenses/>.

MODULE constants  

! Fundamental constants
INTEGER, PARAMETER :: num = KIND(1.D0)
REAL(num),PARAMETER :: pi = 3.1415926535_num
REAL(num),PARAMETER :: mu0 = 4.0D-7*pi

! Background constants
REAL(num),PARAMETER :: rho0 = 1.66D-12
REAL(num),PARAMETER :: B0 = 1.0D-3

! Domain size, time settings and save directory
INTEGER, PARAMETER :: nr=625, nz=2500
REAL(num),PARAMETER :: rmin=0.0_num, rmax=50.0D6
REAL(num),PARAMETER :: zmin=0.0_num, zmax=125.0D6
REAL(num),PARAMETER :: t_end=4500.0_num, t_interval=50.0_num
REAL(num),PARAMETER :: t0=150.0_num !rampup time
CHARACTER(len=100),PARAMETER :: save_dir='Data'

! Setup parameters
REAL(num),PARAMETER :: H = 20.0D6
REAL(num),PARAMETER :: visc = 5.0D7
REAL(num),PARAMETER :: period = 10.0_num
REAL(num),PARAMETER :: alpha = 0.2_num 
REAL(num),PARAMETER :: zeta = 5.0_num
REAL(num),PARAMETER :: u0 = 1.0D5
REAL(num),PARAMETER :: r0 = 5.0D6
REAL(num),PARAMETER :: omega = 2.0_num*pi/period

! Damping
LOGICAL, PARAMETER :: topdamp = .TRUE.
LOGICAL, PARAMETER :: outdamp = .TRUE.

! Restart parameters
! NB If using restart make sure all other parameters are the same as in the previous run.
LOGICAL, PARAMETER :: restart = .TRUE.
CHARACTER(len=100),PARAMETER :: v_in = trim(save_dir)//'/v_060.dat'
CHARACTER(len=100),PARAMETER :: b_in = trim(save_dir)//'/b_060.dat'
REAL(num),PARAMETER :: restart_time = 3000.0_num, last_output = 60
END MODULE constants

MODULE functions  
USE constants
IMPLICIT NONE
CONTAINS

PURE FUNCTION funcv(t,dr,dz,r,z,Br,Bz,rho,v,b)  
! Calculates the RHS of the momentum equation.
REAL(num), INTENT(IN) :: t,dr,dz
REAL(num), dimension (-2:nr+2),INTENT(IN) :: r
REAL(num), dimension (-2:nz+2),INTENT(IN) :: z
REAL(num), dimension (-2:nr+2,-2:nz+2),INTENT(IN) :: Br,Bz,rho,v,b
REAL(num), dimension (-2:nr+2,-2:nz+2) :: dbdr,dbdz,dvdr,dvdz,argdr,argdz,diffdr,diffdz,funcv
INTEGER :: i,j
do i=1,nr
   do j=0,nz
      dbdr(i,j) = (-b(i+2,j)+8.0_num*b(i+1,j)-8.0_num*b(i-1,j)+b(i-2,j))/(12.0_num*dr)
      dbdz(i,j) = (-b(i,j+2)+8.0_num*b(i,j+1)-8.0_num*b(i,j-1)+b(i,j-2))/(12.0_num*dz)
      dvdr(i,j) = (-v(i+2,j)+8.0_num*v(i+1,j)-8.0_num*v(i-1,j)+v(i-2,j))/(12.0_num*dr)
      dvdz(i,j) = (-v(i,j+2)+8.0_num*v(i,j+1)-8.0_num*v(i,j-1)+v(i,j-2))/(12.0_num*dz)
      argdr(i,j) = rho(i,j)*visc*r(i)*dvdr(i,j)
      argdz(i,j) = rho(i,j)*visc*dvdz(i,j)
      diffdr(i,j) = (-argdr(i+2,j)+8.0_num*argdr(i+1,j)-8.0_num*argdr(i-1,j)+argdr(i-2,j))/(12.0_num*dr)
      diffdz(i,j) = (-argdz(i,j+2)+8.0_num*argdz(i,j+1)-8.0_num*argdz(i,j-1)+argdz(i,j-2))/(12.0_num*dz)
      funcv(i,j) = Br(i,j)*dbdr(i,j)+Bz(i,j)*dbdz(i,j)+(Br(i,j)*b(i,j)/r(i))
      funcv(i,j) = (1.0_num/(mu0*rho(i,j)))*funcv(i,j)
      funcv(i,j) = funcv(i,j)+(1.0_num/(r(i)*rho(i,j)))*diffdr(i,j)+(1.0_num/rho(i,j))*diffdz(i,j)
   enddo
enddo
END FUNCTION funcv

PURE FUNCTION funcb(t,dr,dz,r,z,Br,Bz,rho,v,b)   
! Calculates the RHS of the induction equation.
REAL(num), INTENT(IN) :: t,dr,dz
REAL(num), dimension (-2:nr+2),INTENT(IN) :: r
REAL(num), dimension (-2:nz+2),INTENT(IN) :: z
REAL(num), dimension (-2:nr+2,-2:nz+2),INTENT(IN) :: Br,Bz,rho,v,b
REAL(num), dimension (-2:nr+2,-2:nz+2) :: dvdr,dvdz,funcb
INTEGER :: i,j
do i=1,nr
   do j=0,nz
      dvdr(i,j) = (-v(i+2,j)+8.0_num*v(i+1,j)-8.0_num*v(i-1,j)+v(i-2,j))/(12.0_num*dr)
      dvdz(i,j) = (-v(i,j+2)+8.0_num*v(i,j+1)-8.0_num*v(i,j-1)+v(i,j-2))/(12.0_num*dz)
      funcb(i,j) = (Br(i,j)*dvdr(i,j)+Bz(i,j)*dvdz(i,j)-(Br(i,j)*v(i,j)/r(i)))
   enddo
enddo
END FUNCTION funcb

SUBROUTINE edgebc(f) 
! Boundary conditions for domain edges 
!(i,e. for radial boundaries)
USE constants
REAL(num), dimension (-2:nr+2,-2:nz+2) :: f
! radial antisymmetry for the r=0 boundary makes sense, 
! for the outer radial boundary we use reflective BC.
f(-2,:)=-f(2,:) 
f(-1,:)=-f(1,:)      
f(nr+1,:)=f(nr-1,:)
f(nr+2,:)=f(nr-2,:)
END SUBROUTINE edgebc

SUBROUTINE drivebc_v(f,r,z,Va,t)  
! Velocity BC for lower boundary driving
USE constants
REAL(num), dimension (-2:nr+2,-2:nz+2) :: f,Va
REAL(num), dimension (-2:nr+2) :: r
REAL(num), dimension (-2:nz+2) :: z
REAL(num) :: t
INTEGER :: i,j
do i=0,nr
   do j=-2,-1
      if (r(i) .LT. r0) then
         f(i,j)=u0*(r(i)/r0)*(1.0_num-(r(i)/r0)**2) & !Amplitude
              *sin(omega*(z(j)/Va(i,j)-t)) & !Phase
              *(1.0_num-exp(-(t/t0)**3)) !Rampup
      else
         f(i,j) = 0.0_num
      endif
   enddo
enddo
END SUBROUTINE drivebc_v

SUBROUTINE drivebc_b(f,r,z,Va,rho,t)  
! Magnetic BC for lower boundary driving
USE constants
 REAL(num), dimension (-2:nr+2,-2:nz+2) :: f,Va,rho
REAL(num), dimension (-2:nr+2) :: r
REAL(num), dimension (-2:nz+2) :: z
REAL(num) :: t
INTEGER :: i,j
do i=0,nr
   do j=-2,-1
      if (r(i) .LT. r0) then
         f(i,j)=-u0*(r(i)/r0)*(1.0_num-(r(i)/r0)**2)*sqrt(mu0*rho(i,j)) & !Amplitude
              *sin(omega*(z(j)/Va(i,j)-t)) & !Phase
              *(1.0_num-exp(-(t/t0)**3)) !Rampup
      else
         f(i,j) = 0.0_num
      endif
   enddo
enddo
END SUBROUTINE drivebc_b

SUBROUTINE damp_top(f,z,dt)  
! Damping BC for upper boundary
USE constants
REAL(num), dimension (-2:nr+2,-2:nz+2) :: f
REAL(num), dimension (-2:nz+2) :: z
REAL(num) :: d,a,dt
INTEGER :: i,j
d = 4.0_num*zmax/5.0_num
do i=-2,nr+2
   do j=-2,nz+2
      if (z(j) .GT. d) then
         a = 1.0_num+dt*(z(j)-d)/(zmax-d)
         f(i,j)= f(i,j)/a
      endif
   enddo
enddo
END SUBROUTINE damp_top

SUBROUTINE dampout(f,r,dt)  
! Damping BC for upper boundary
USE constants
REAL(num), dimension (-2:nr+2,-2:nz+2) :: f
REAL(num), dimension (-2:nr+2) :: r
REAL(num) :: d,a,dt
INTEGER :: i,j
d = 4.0_num*rmax/5.0_num
do i=-2,nr+2
   do j=-2,nz+2
      if (r(i) .GT. d) then
         a = 1.0_num+dt*(r(i)-d)/(rmax-d)
         f(i,j)= f(i,j)/a
      endif
   enddo
enddo
END SUBROUTINE dampout


SUBROUTINE envelope(env_v,env_b,v,b)  
! Update envelope variables
USE constants
REAL(num), dimension (-2:nr+2,-2:nz+2) :: env_v,env_b,v,b
INTEGER :: i,j
do i=-2,nr+2
   do j=-2,nz+2
      env_v(i,j) = max(env_v(i,j),abs(v(i,j)))
      env_b(i,j) = max(env_b(i,j),abs(b(i,j)))
   enddo
enddo
END SUBROUTINE envelope

END MODULE functions

PROGRAM main
USE constants
USE functions
IMPLICIT NONE

REAL(num) :: time,out_time,runtime
REAL(num) :: dr,dz,dt,psib
REAL(num), dimension (-2:nr+2) :: r
REAL(num), dimension (-2:nz+2) :: z
REAL(num), dimension (-2:nr+2,-2:nz+2) :: Br,Bz,phi,psi,rho,Va,v,b,env_v,env_b
REAL(num), dimension (-2:nr+2,-2:nz+2) :: fk1,fm1,fk2,fm2,fk3,fm3,k4,m4
INTEGER :: i,j,step,outname
CHARACTER(len=14) :: filename

! Output Header
print*,""
print*,"    _      _   __    ______     ______     __        _____    _      _    _____    _     _    _____  " 
print*,"   /_/\  /\_\ /\_\  /_/\___\   /_/\___\   /\_\     /\_____\  /_/\  /\_\  /\___/\  /_/\ /\_\ /\_____\ "
print*,"   ) ) )( ( ( \/_/  ) ) ___/   ) ) ___/  ( ( (    ( (_____/  ) ) )( ( ( / / _ \ \ ) ) ) ( (( (_____/ "
print*,"  /_/ //\\ \_\ /\_\/_/ /  ___ /_/ /  ___  \ \_\    \ \__\   /_/ //\\ \_\\ \(_)/ //_/ / \ \_\\ \__\   "
print*,"  \ \ /  \ / // / /\ \ \_/\__\\ \ \_/\__\ / / /__  / /__/_  \ \ /  \ / // / _ \ \\ \ \_/ / // /__/_  "
print*,"   )_) /\ (_(( (_(  )_)  \/ _/ )_)  \/ _/( (_____(( (_____\  )_) /\ (_(( (_( )_) )\ \   / /( (_____\ "
print*,"   \_\/  \/_/ \/_/  \_\____/   \_\____/   \/_____/ \/_____/  \_\/  \/_/ \/_/ \_\/  \_\_/_/  \/_____/ "
print*,""                                                                                                    

! Calculate grid spacings and tube boundary
dr=(rmax-rmin)/real(nr)
dz=(zmax-zmin)/real(nz)
psib = (r0**2)/(2*H)

! Calculate background values over grid
! These are r,z,Br,Bz,psi,rho and Va
do i=-2,nr+2
   do j=-2,nz+2
      r(i) = abs(i*dr)
      z(j) = j*dz
      Br(i,j) = B0*exp(-z(j)/H)*bessel_j1(r(i)/H)
      Bz(i,j) = B0*exp(-z(j)/H)*bessel_j0(r(i)/H)
      phi(i,j) = -H*exp(-z(j)/H)*bessel_j0(r(i)/H)
      psi(i,j) = r(i)*exp(-z(j)/H)*bessel_j1(r(i)/H)
      if (psi(i,j) .LE. psib) then
         rho(i,j) = (rho0/zeta)*(1.0_num+(zeta-1.0_num)*((1.0_num-(psi(i,j)/psib))**2))
      else
         rho(i,j) = (rho0/zeta)
      endif
      rho(i,j) = rho(i,j)*exp(-alpha*z(j)/H)
      Va(i,j) = sqrt((Br(i,j)**2)+(Bz(i,j)**2))/sqrt(mu0*rho(i,j))
   enddo
enddo

! Save background values
OPEN(unit=60, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/Va.dat')
write(60) Va(-2:nr+2,-2:nz+2)
close(60)
OPEN(unit=61, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/Br.dat')
write(61) Br(-2:nr+2,-2:nz+2)
close(61)
OPEN(unit=62, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/Bz.dat')
write(62) Bz(-2:nr+2,-2:nz+2)
close(62)
OPEN(unit=63, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/rho.dat')
write(63) rho(-2:nr+2,-2:nz+2)
close(63)
OPEN(unit=64, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/phi.dat')
write(64) phi(-2:nr+2,-2:nz+2)
close(64)
OPEN(unit=65, FORM = 'UNFORMATTED', FILE = trim(save_dir)//'/psi.dat')
write(65) psi(-2:nr+2,-2:nz+2)
close(65)

! Calculate timestep
! This is based on the von neumann stability conditions for our PDEs
dt=0.1_num*min(min(dr,dz)/maxval(Va),min(dr,dz)**2/visc)

print '("dr = ",F17.0)', dr
print '("dz = ",F17.0)', dz
print '("dt = ",F8.4," seconds.")',dt

! Initialise time
if (restart) then
   step = 0
   out_time = restart_time + t_interval
   time = restart_time
   outname = last_output + 1
else
   step = 0
   out_time = t_interval
   time = 0.0_num
   outname = 1
endif

! Initialise v and b
if (restart) then
   open(unit=1,form='unformatted',status='old',action='read',file=v_in)
   read(1) v
   close(1)
   open(unit=2,form='unformatted',status='old',action='read',file=b_in)
   read(2) b
   close(2)
else
   v(:,:) = 0.0_num
   b(:,:) = 0.0_num
endif

! Initialise RK4 parameters
fk1(:,:) = 0.0_num
fm1(:,:) = 0.0_num
fk2(:,:) = 0.0_num
fm2(:,:) = 0.0_num
fk3(:,:) = 0.0_num
fm3(:,:) = 0.0_num
k4(:,:) = 0.0_num
m4(:,:) = 0.0_num

! Main loop

do while (time .le. t_end) 

! RK4 for a system of two coupled equations
!
! k1,k2,k3,m1,m2,m3 are not calculated directly
! fk1 = v+k1/2     fm1 = b+m1/2
! fk2 = v+k2/2     fm2 = b+m2/2
! fk3 = v+k3       fm3 = b+m3
! This is because the BC must be applied to functions and not differences.
! BC must be applied to fk1,fk2,fk3,fm1,fm2,fm3, this is because derivatives of them are taken.
! No BC are required for k4 and m4 because no derivatives of them are taken.

fk1 = v+(dt*funcv(time,dr,dz,r,z,Br,Bz,rho,v,b))/2.0_num
fm1 = b+(dt*funcb(time,dr,dz,r,z,Br,Bz,rho,v,b) )/2.0_num
call drivebc_v(fk1,r,z,Va,time+dt/2.0_num)
call drivebc_b(fm1,r,z,Va,rho,time+dt/2.0_num)
call edgebc(fk1)
call edgebc(fm1)

fk2 = v+(dt*funcv(time+dt/2.0_num,dr,dz,r,z,Br,Bz,rho,fk1,fm1))/2.0_num 
fm2 = b+(dt*funcb(time+dt/2.0_num,dr,dz,r,z,Br,Bz,rho,fk1,fm1))/2.0_num
call drivebc_v(fk2,r,z,Va,time+dt/2.0_num)
call drivebc_b(fm2,r,z,Va,rho,time+dt/2.0_num)
call edgebc(fk2)
call edgebc(fm2)

fk3 = v+dt*funcv(time+dt/2.0_num,dr,dz,r,z,Br,Bz,rho,fk2,fm2)
fm3 = b+dt*funcb(time+dt/2.0_num,dr,dz,r,z,Br,Bz,rho,fk2,fm2)
call drivebc_v(fk3,r,z,Va,time+dt)
call drivebc_b(fm3,r,z,Va,rho,time+dt)
call edgebc(fk3)
call edgebc(fm3)

k4 = dt*funcv(time,dr,dz,r,z,Br,Bz,rho,fk3,fm3) 
m4 = dt*funcb(time,dr,dz,r,z,Br,Bz,rho,fk3,fm3)

v = v + (2.0_num*(fk1-v)+4.0_num*(fk2-v)+2.0_num*(fk3-v)+k4)/6.0_num
b = b + (2.0_num*(fm1-b)+4.0_num*(fm2-b)+2.0_num*(fm3-b)+m4)/6.0_num
call drivebc_v(v,r,z,Va,time+dt)
call drivebc_b(b,r,z,Va,rho,time+dt)
call edgebc(v)
call edgebc(b)
if (topdamp) then
   call damp_top(v,z,dt)
   call damp_top(b,z,dt)
endif
if (outdamp) then
   call dampout(v,r,dt)
   call dampout(b,r,dt)
endif

call envelope(env_v,env_b,v,b)

! Increment timestep
time =time + dt
step = step + 1

! Output block
if (time-dt.LT.out_time .AND. time.GT.out_time) then
print*,''
! Save data for v and b
write (filename, "(A3,I0.3,A4)") "/v_", outname,'.dat'
OPEN(unit=66, FORM = 'UNFORMATTED', FILE = trim(save_dir)//filename)
write(66) v(-2:nr+2,-2:nz+2)
close(66)
write (filename, "(A3,I0.3,A4)") "/b_", outname,'.dat'
OPEN(unit=67, FORM = 'UNFORMATTED', FILE = trim(save_dir)//filename)
write(67) b(-2:nr+2,-2:nz+2)
close(67)
! Save envelope data
write (filename, "(A7,I0.3,A4)") "/env_v_", outname,'.dat'
OPEN(unit=68, FORM = 'UNFORMATTED', FILE = trim(save_dir)//filename)
write(68) env_v(-2:nr+2,-2:nz+2)
close(68)
write (filename, "(A7,I0.3,A4)") "/env_b_", outname,'.dat'
OPEN(unit=69, FORM = 'UNFORMATTED', FILE = trim(save_dir)//filename)
write(69) env_b(-2:nr+2,-2:nz+2)
close(69)
! Readouts for simulation time, step and runtime
call cpu_time(runtime)
print '("output",I3,".")', outname
print '("simulation time = ",F8.2,".")', time
print '("step number = ",I12,".")', step
print '("runtime = ",F8.1," seconds.")',runtime
! Increment counters
out_time = out_time + t_interval
outname = outname + 1
endif

enddo

END PROGRAM main

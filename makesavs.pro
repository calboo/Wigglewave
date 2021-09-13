pro makesavs

nr = 500
nz = 2500
dir = 'Data'
snapshot = '060'

va = dblarr(nr+5,nz+5) 
OPENR, 1, dir+'/Va.dat',  /F77_UNFORMATTED
READU, 1, va
CLOSE, 1 

Br = dblarr(nr+5,nz+5)  
OPENR, 2, dir+'/Br.dat',  /F77_UNFORMATTED
READU, 2, Br
CLOSE, 2

Bz = dblarr(nr+5,nz+5)  
OPENR, 3, dir+'/Bz.dat',  /F77_UNFORMATTED
READU, 3, Bz
CLOSE, 3

rho = dblarr(nr+5,nz+5)  
OPENR, 4, dir+'/rho.dat',  /F77_UNFORMATTED
READU, 4, rho
CLOSE, 4

phi = dblarr(nr+5,nz+5) 
OPENR, 5, dir+'/phi.dat',  /F77_UNFORMATTED
READU, 5, phi
CLOSE, 5

psi = dblarr(nr+5,nz+5) 
OPENR, 6, dir+'/psi.dat',  /F77_UNFORMATTED
READU, 6, psi
close, 6

v = dblarr(nr+5,nz+5)  
OPENR, 7, dir+'/v_'+snapshot+'.dat',  /F77_UNFORMATTED
READU, 7, v
CLOSE, 7

b = dblarr(nr+5,nz+5)  
OPENR, 8, dir+'/b_'+snapshot+'.dat',  /F77_UNFORMATTED
READU, 8, b
CLOSE, 8

env_v = dblarr(nr+5,nz+5)  
OPENR, 9, dir+'/env_v_'+snapshot+'.dat',  /F77_UNFORMATTED
READU, 9, env_v
CLOSE, 9

env_b = dblarr(nr+5,nz+5)  
OPENR, 10, dir+'/env_b_'+snapshot+'.dat',  /F77_UNFORMATTED
READU, 10, env_b
CLOSE, 10

save, va, filename= dir+'/Va.sav'
save, Br, filename= dir+'/Br.sav'
save, Bz, filename= dir+'/Bz.sav'
save, rho, filename= dir+'/rho.sav'
save, phi, filename= dir+'/phi.sav'
save, psi, filename= dir+'/psi.sav'
save, v, filename= dir+'/v.sav'
save, b, filename= dir+'/b.sav'
save, env_v, filename= dir+'/env_v.sav'
save, env_b, filename= dir+'/env_b.sav'

end

pro wave_energy

;;  Define parameters
;;  NB don't include damping region
mu0 = 4.0d0*!dpi*1.0d-7
H = 20.0d6
b0 = 1.0d-3
rho0 = 1.0d0*1.66d-12
alpha = 0.2d0
rsize = 500
zsize = 2000 
rdim = 40.0d6
zdim = 100.0d6
zeta = 5.0d0
u0 = 1.0d3
r0 = 5.0d6

save_dir= 'Data'
run_name= ''
save_loc = save_dir+run_name

;;  Initialise arrays
r = findgen(rsize+1)*rdim/rsize
z = findgen(zsize+1)*zdim/zsize
waven = dblarr(rsize+1,zsize+1)

;; Resore data
restore, save_loc+'/rho.sav'
restore, save_loc+'/phi.sav'
restore, save_loc+'/psi.sav'
restore, save_loc+'/env_v.sav'
restore, save_loc+'/env_b.sav'

rho = rho[2:rsize+2,2:zsize+2]
phi = phi[2:rsize+2,2:zsize+2]
psi = psi[2:rsize+2,2:zsize+2]
env_v = env_v[2:rsize+2,2:zsize+2]
env_b = env_b[2:rsize+2,2:zsize+2]

;;  Calculate psi at boundary
psib = (r0^2.0)/(2.0*H)

;; Calculate integrand for wave energy
for i = 0,rsize do begin
   for j = 0,zsize do begin
      waven(i,j) = ((!dpi*b0*H)/mu0)*(env_v[i,j]*env_b[i,j])
   endfor
endfor

;window, 0
;shade_surf, waven
save, waven, filename= save_loc+'/waven.sav'

;; Plot Integrand for wave energy and field lines
window, 1, xsize=800, ysize=600
contour,  alog(waven+1.0d-12), /fill, nlevels=100

nlvl=40
lvlphi= fltarr(nlvl)
for i= 0,nlvl-1 do begin
   lvlphi(i) = min(phi)+(i^0.5d0)*(max(phi)-min(phi))/((nlvl-1)^0.5d0)
endfor
lvlpsi= fltarr(nlvl)
for i= 0,nlvl-1 do begin
   lvlpsi(i) = min(psi)+(i^2.0d0)*(max(psi)-min(psi))/((nlvl-1)^2.0d0)
endfor
contour, psi,levels=lvlpsi, /overplot
contour, phi,levels=lvlphi, /overplot  

;;  Integrate wave energy across field lines at each height 
;;  (also overplots paths for each integration)

zscale = [] 
hscale = []
en_lvl = []
z0 = -H*alog(beselj(r0/H,0))
zlevels = zsize

for k= 0,zlevels-1 do begin
z_en = k*zdim/zlevels
if (z_en LT z0) then begin 
   print, string(100*float(k)/float(zlevels),FORMAT='(f4.1)')+' %',' too low'
endif else begin
   phi1 = -H*exp(-z_en/H)
   i = 0
   pointi =[]
   pointj =[]
   psi_int = []
   waven_int = []
   while (i LT rsize) do begin
      p1 = value_locate(phi[i,*],phi1)
      p2 = p1+1
      if (p2 GT zsize) then break
      psi_i = psi(i,p1)+(phi1-phi(i,p1))*(psi(i,p2)-psi(i,p1))/(phi(i,p2)-phi(i,p1))
      waven_i = waven(i,p1)+(phi1-phi(i,p1))*(waven(i,p2)-waven(i,p1))/(phi(i,p2)-phi(i,p1))
      if (psi_i GE psib) then break
      pointi = [pointi, i]
      pointj = [pointj,p2]
      waven_int = [waven_int, waven_i]
      psi_int = [psi_int, psi_i]
      i = i+1 
   endwhile
   print, string(100*float(k)/float(zlevels),FORMAT='(f4.1)')+' %'
   oplot, pointi,pointj
   en1 = int_tabulated(psi_int,waven_int)
   zscale = [zscale, z_en]
   hscale = [hscale, z_en/H]
   en_lvl = [en_lvl, en1]
   endelse
endfor 
en_lvl_norm = en_lvl/(en_lvl[0])

; Plot total wave energy 

p1 = plot(zscale/1e6, en_lvl/1e15, THICK=2, LINESTYLE=0, $
          XTITLE= 'Height / Mm', YTITLE= 'Wave Energy Flux / PW', $
         FONT_SIZE=15);, YRANGE=[0.0,2.7])
p2 = plot(zscale/1e6, en_lvl_norm, THICK=2, LINESTYLE=0, $
          XTITLE= 'Height / Mm', YTITLE= 'Normalised Wave Energy Flux', $
         FONT_SIZE=15, YRANGE=[0.0,1.0])

save, zscale, filename= save_loc+'/zscale.sav'
save, hscale, filename= save_loc+'/hscale.sav'
save, en_lvl, filename= save_loc+'/en_lvl_0.sav'
save, en_lvl_norm, filename=save_loc+'/en_lvl_norm_0.sav'

end


pro plotgraphs

;; DESCRIPTION

; IDL script used to produce Figures 10 and 11 in the paper,
; Enhanced phase mixing of torsional Alfven waves in stratified and
; divergent solar coronal structures – Paper I. Linear solutions
; C.Boocock and D.Tsiklauri

; Two panel plots are produced using the outputs from TAWAS
; and Wigglewave. The first plot compares the azimuthal velocity 
; between the two outputs and the second plot compares the 
; azimuthal magnetic field perturbation between the two outputs.

; The first panel shows a plot of the azimuthal velocity from the
; TAWAS output, the azimuthal velocity from the Wigglewave output
; and the difference between the wave envelopes for the azimuthal
; velocity between TAWAS and Wigglewave.

; The second panel shows a plot of the azimuthal magnetic field
; perturbation from the TAWAS output, the azimuthal magnetic field
; perturbation from the Wigglewave output and the difference between
; the wave envelopes for the azimuthal magnetic field perturbation
; between TAWAS and Wigglewave.

;; USAGE

; This script must be run in IDL, to run this script correctly
; outputs from TAWAS and Wigglewave that use the same
; parameters should be compared. To make this simple the outputs
; should be saved under identical path structures under seperate
; root directories for TAWAS and Wigglewave.

; The contour levels will need to be adjusted depending
; on the parameters being considered.

; Figures 10 and 11 were produced by using this script
; with outputs generated from the following parameters:
; H = 50 Mm, H_rho= 50 Mm, T = 10 s,
;u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.

;; SCRIPT

; Specify directories for input

dir = 'WKB_study/Data/T60/'
anadir = '../Analytic/'+dir
wigdir = '../Wigglewave/'+dir

; Set up the coordinates

r = 15*findgen(501)/500
z = 100*findgen(2001)/2000

; Azimuthal velocity comparison panel

; Restore data from TAWAS
restore, anadir+'v.sav'
v = real_part(v)
print,'velocity range TAWAS'
print, max(v), min(v)
print, ''

; Plot the first contour plot using data from TAWAS
ct = colortable(5)
clevels = 60*(findgen(101)/50-1.0)
w1 = WINDOW(DIMENSIONS=[1800,600])
c1 = CONTOUR(v/1e3,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',YTITLE='z / Mm',FONT_SIZE=20,$
             C_VALUE=clevels,POSITION=[0.07,0.25,0.35,0.95],/CURRENT)
cb = COLORBAR(TARGET=c1,TITLE='Velocity / kms!E-1',$
              POSITION=[0.05,0.12,0.35,0.14],FONT_SIZE=20)

; Restore data from Wigglewave
restore, wigdir+'v.sav'
v = v[2:502,2:2002]
v = real_part(v)
print,'velocity range Wigglewave'
print, max(v), min(v)
print, ''

; Plot the second contour plot using data from Wigglewave
ct = colortable(5)
clevels = 60*(findgen(101)/50-1.0)
c2 = CONTOUR(v/1e3,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',FONT_SIZE=20,$
             C_VALUE=clevels,POSITION=[0.39,0.25,0.67,0.95],/CURRENT)
cb = COLORBAR(TARGET=c2,TITLE='Velocity / kms!E-1',$
              POSITION=[0.37,0.12,0.67,0.14],FONT_SIZE=20)

; Restore the envelope data for the azimuthal velocity
; from both TAWAS and Wigglewave
restore, anadir+'v_env.sav'
restore, wigdir+'env_v.sav'
env_v = env_v[2:502,2:2002]

; Calculate the difference between the two envelopes
v_diff = abs(env_v-v_env)
print,'velocity range difference'
print, max(v_diff), min(v_diff) 
print, ''

; Plot the third contour plot showing the difference 
ct = colortable(3)
clevels = 2*(findgen(101)/100)
c3 = CONTOUR(v_diff/1e3,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',FONT_SIZE=20,$
             C_VALUE=clevels,POSITION=[0.71,0.25,0.99,0.95],/CURRENT)
cb = COLORBAR(TARGET=c3,TITLE='Difference / kms!E-1',$
              POSITION=[0.69,0.12,0.99,0.14],FONT_SIZE=20)

; Azimuthal magnetic field comparison panel

; Restore data from TAWAS
restore, anadir+'b.sav'
b = real_part(bwave)
print,'magnetic field range TAWAS'
print, max(b), min(b)
print, ''

; Plot the first contour plot using data from TAWAS
ct = colortable(5)
clevels = 5e1*(findgen(101)/50-1.0)
w4 = WINDOW(DIMENSIONS=[1800,600])
c4 = CONTOUR(b*1e6,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',YTITLE='z / Mm',FONT_SIZE=20,$
             C_VALUE=clevels,POSITION=[0.07,0.25,0.35,0.95], /CURRENT)
cb = COLORBAR(TARGET=c4,TITLE='Magnetic field / $\mu$T',$
              POSITION=[0.05,0.12,0.35,0.14],FONT_SIZE=20)

; Restore data from Wigglewave
restore, wigdir+'b.sav'
b = b[2:502,2:2002]
b = real_part(b)
print,'magnetic field range Wigglewave'
print, max(b), min(b)
print, ''

; Plot the second contour plot using data from Wigglewave
ct = colortable(5)
clevels = 5e1*(findgen(101)/50-1.0)
c5 = CONTOUR(b*1e6,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',FONT_SIZE=20,$
              C_VALUE=clevels,POSITION=[0.39,0.25,0.67,0.95],/CURRENT)
cb = COLORBAR(TARGET=c5,TITLE='Magnetic field / $\mu$T',$
              POSITION=[0.37,0.12,0.67,0.14],FONT_SIZE=20)

; Restore the envelope data for the azimuthal magnetic field
; from both TAWAS and Wigglewave
restore, anadir+'b_env.sav'
restore, wigdir+'env_b.sav'
env_b = env_b[2:502,2:2002]

; Calculate the difference between the two envelopes
b_diff = abs(env_b-b_env)
print,'magnetic field range difference'
print, max(b_diff), min(b_diff) 
print, ''

; Plot the third contour plot showing the difference 
ct = colortable(3)
clevels = (findgen(101)/100)
c6 = CONTOUR(b_diff*1e6,r,z, SHADING=1, /FILL, RGB_TABLE=ct, N_LEVELS=100,$
             XTITLE='r / Mm',FONT_SIZE=20,$
             C_VALUE=clevels,POSITION=[0.71,0.25,0.99,0.95],/CURRENT)
cb = COLORBAR(TARGET=c6,TITLE='Difference / $\mu$T',$
              POSITION=[0.69,0.12,0.99,0.14],FONT_SIZE=20)

; Save both panels as png files

c1.save, "Images/v_panel.png"
c4.Save, "Images/b_panel.png"

end

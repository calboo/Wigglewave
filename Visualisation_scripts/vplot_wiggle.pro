pro vplot_wiggle

;; DESCRIPTION

; IDL script used to produce Figure 3 in the paper,
; Enhanced phase mixing of torsional Alfven waves in stratified and
; divergent solar coronal structures – Paper II. Non-linear simulations
; C.Boocock and D.Tsiklauri

; Produces a graph of the azimuthal velocity against the radius, r,
; and height, z, for the Wigglewave output specified in the script.

;; USAGE

; This script must be run in IDL, the user must specify a valid filepath
; within the script corresponding to a directory containing Wigglewave
; outputs.

; The contour levels will need to be adjusted depending
; on the parameters being considered.

; Figure 3 was produced by using this script with Wigglewave
; outputs generated from the following parameters:
; H = 50 Mm, H_rho= 50 Mm, T = 60 s,
; u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.

;; SCRIPT

; Specify directory for input

dir = 'WKB_study/Data/T60/'

; Read azimuthal velocity data from Wigglewave outputs

restore, dir+'v.sav'
v = v[2:502,2:2002]
v = real_part(v)
print, max(v), min(v)

; Set up coordinates for contour plot

r = 15*findgen(501)/500
z = 100*findgen(2001)/2000

; Plot a coloured contour of the azimuthal velocity

w1 = WINDOW(DIMENSIONS=[800,800])
ct = colortable(3)
clevels = 100.0*(findgen(101)/50-1.0)
c1 = CONTOUR(v/1e3,r,z, SHADING=1, /FILL, RGB_TABLE=ct,$
             XTITLE='r / Mm',YTITLE='z / Mm', FONT_SIZE=20,$
             XSTYLE=1,YSTYLE=1,/CURRENT,$
             C_VALUE=clevels,POSITION=[0.12,0.25,0.95,0.95])
cb = COLORBAR(TARGET=cw,TITLE='Velocity / kms!E-1',$
              POSITION=[0.05,0.14,0.95,0.16],FONT_SIZE=20)

; Save image as a png file

c1.save, "vplot_wiggle.png"

end

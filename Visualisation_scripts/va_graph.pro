pro va_graph

;; DESCRIPTION

; IDL script used to produce Figures 1 and 2 in the paper,
; Enhanced phase mixing of torsional Alfven waves in stratified and
; divergent solar coronal structures – Paper II. Non-linear simulations
; C.Boocock and D.Tsiklauri

; Two plots are produced using the specified output data from Wigglewave. 

; The first plot is a graph showing how the the normalised
; values for the equilibrium density and Alfven speed vary against
; radius across the lower boundary of the domain.

; The second plot is a coloured contour plot of the Alfven speed
; plotted against the radius, r, and height, z. This is effectively
; a contour of the Alfven speed across the whole domain due to the 
; condition of axisymmetry in the initial conditions.

;; USAGE

; This script must be run in IDL, the user must specify valid filepaths
; within the script corresponding to Wigglewave outputs for the Alfven
; speed and density.

; The contour levels will need to be adjusted depending
; on the parameters being considered.

; Figures 1 and 2 were produced by using this script with Wigglewave
; outputs generated from the following parameters:
; H = 50 Mm, H_rho= 50 Mm, T = 60 s,
; u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.

;; SCRIPT

; Specify directories for input

restore, 'WKB_study/Data/T60/Va.sav'
restore, 'WKB_study/Data/T60/rho.sav'

; Set up the coordinates

r = findgen(501)/501*15.03
z = findgen(2501)/2501*100.04
va = va[2:502,2:2502]

; Set up character of accented e
eaccent = string("351B)
;; This line is to correct code colouring "

; Plot a graph of the normalised Alfven speed
; at the lower boundary against the radius.

p1 = plot(r,va[*,0]/692374,LINESTYLE=2,YRANGE=[0,2.5],THICK=2, $
         YTITLE = 'Normalised $\rho$!D0!N and V!DA!N', $
         XTITLE = 'r / Mm')

; Overplot a graph of the normalised density
; at the lower boundary against the radius

p2 = plot(r,rho[*,0]/1.66e-12,LINESTYLE=0,THICK=2,/OVERPLOT)

; Save this first plot as a png file 

p1.Save, "rho_va.png"
 
; Set up a seperate window 

w1 = WINDOW(DIMENSIONS=[800,800])

; Plot a coloured contour of the equilibrium Alfven speed
; against the radius, r, and height, z

ct = colortable(3)
clevels = 1000.0*(findgen(101)/100)
c1 = CONTOUR(va/1e3,r,z, SHADING=1,$
             /FILL, RGB_TABLE=ct, XTITLE='r / Mm',YTITLE='z / Mm',$
             XSTYLE=1,YSTYLE=1,N_LEVELS=100,/CURRENT,$
             POSITION=[0.12,0.25,0.95,0.95],FONT_SIZE=20,C_VALUE=clevels)
cb = COLORBAR(TARGET=c1,TITLE='Alfv'+eaccent+'n Speed / kms!U-1',$
              POSITION=[0.05,0.14,0.95,0.16],FONT_SIZE=20)

; Save this second plot to a png file

c1.Save, "va.png"

end

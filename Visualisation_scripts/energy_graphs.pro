pro energy_graphs

;; DESCRIPTION

; IDL script used to produce Figures 12, 13 and 14 in the paper,
; Enhanced phase mixing of torsional Alfven waves in stratified and
; divergent solar coronal structures – Paper I. Linear solutions
; C.Boocock and D.Tsiklauri

; Plots of the normalised wave energy flux across a given magnetic
; surface, against the height at which that surface intersects the 
; vertical axis, are plotted. Each plot has a line for the normalised
; wave energy flux as calculated using TAWAS (a solid line) and a line
; for the normalised wave energy flux as calculated using Wigglewave
; (a dashed line).

;; USAGE

; This script must be run in IDL, to run this script correctly
; the outputs from TAWAS and Wigglewave that use the same
; parameters should be compared. To make this simple the outputs
; should be saved under identical path structures under seperate
; root directories for TAWAS and Wigglewave.

; Figures 12 was produced by using this script with outputs
; generated from the following parameters:
; H = 50 Mm, H_rho= 50 Mm, T = 10 s,
; u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.

; Figures 13 was produced by using this script with outputs
; generated from the following parameters:
; H = 20 Mm, H_rho= 100 Mm, T = 10 s,
; u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.

; Figures 14 was produced by using this script with outputs
; generated from the following parameters:
; H = 50 Mm, H_rho= 50 Mm, T = 60 s,
; u0 = 100 kms-1 and viscosity= 5 ×10^7 m2s−1.
; This figure was then further annotated using powerpoint to show
; the effects of reflection and dissipation at different heights.

;; SCRIPT

; Specify directories for input

dir = 'WKB_study/Data/T60/'
anadir = '../Analytic/'+dir
wigdir = '../Wigglewave/'+dir

; Restore normalised wave energy flux
; and height data from Wigglewave

restore, wigdir+'/en_lvl_norm_0.sav'
restore, wigdir+'/zscale.sav'

; Plot normalised wave energy flux from Wigglewave data

p1 = plot(zscale/1e6, en_lvl_norm, THICK=2, LINESTYLE=2, $
          XTITLE= 'Height / Mm', YTITLE= 'Normalised Wave Energy Flux', $
         FONT_SIZE=15, YRANGE=[0.0,1.0])

; Restore normalised wave energy flux from TAWAS

restore, anadir+'/en_lvl_norm_0.sav'

; Overplot normalised wave energy flux from TAWAS data

p1 = plot(zscale/1e6, en_lvl_norm, THICK=2, LINESTYLE=0, $
          XTITLE= 'Height / Mm', YTITLE= 'Normalised Wave Energy Flux', $
         FONT_SIZE=15,/overplot)

; Save graph image as a png file

p1.save, "Images/power_T60.png"

end

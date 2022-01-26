; analyze excess TiN QP.pro
; 
; generates histograms and x/y-plots of pulses fitted by 'pulse fitting excess TiN QP.pro'
; -->> these pulses are all identified membrane hits (of only one channel) and are fitted asuming a thermal and a non-thermal signal component
; saves histograms & plots in 'sample+date measurement channel 1/2 analyzed (Excess QP).ps'
;
; <><><>  adjust before running:       <><><>
; <><><>  information in 'to_fit.txt'  <><><>
; <><><>  channel to analyze           <><><>
; <><><>  MaxChisq & MinRiseTime       <><><>
; <><><>  remove histogram bounds      <><><>

;----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  channel = 2  ; channel to analyze
  scatter = 0   ; do scatter plots or not
  summary = 1   ; save summatry file or not
  
  MinRiseTime = 150.0
  MaxChisq = 100.0
;<><><><><><>

; masterfile that defines where to find the data to analyze:
path = '/home/gerhard/' ; where to find the file 'to_fit.txt'
cd,path
masterfile = 'to_fit.txt'
line=string(200)
openr,1,strcompress(path+masterfile)
readf,1,line & datapath = (strsplit(line,/EXTRACT))[0] ; read path where to find the data file
readf,1,line & outpath = (strsplit(line,/EXTRACT))[0]  ; read path where to save output file
readf,1,line & pulsename = (strsplit(line,/EXTRACT))[0]          ; file with pulse data to be analyzed
readf,1,line & StartRes = (strsplit(line,'resnum=',/EXTRACT))[0] ; resonator group as numbered by the IQ sweep routine
readf,1,line & atten = (strsplit(line,'atten=',/EXTRACT))[0]     ; unused: attenuator setting for pulse file in dB
readf,1,line & temperature = (strsplit(line,'T=',/EXTRACT))[0]   ; measurement temperature in mK
readf,1,line & temp = (strsplit(line,'sub=',/EXTRACT))[0]        ; unused: do a background subtraction or not
readf,1,line & temp = (strsplit(line,'loss=',/EXTRACT))[0]       ; unused: attenuation before device, DilutionFridge = 45, ADR = 35
readf,1,line & Npts = (strsplit(line,'Np=',/EXTRACT))[0]         ; number of points in each pulse
readf,1,line & sample = (strsplit(line,/EXTRACT))[0]             ; read sample name
close,1

; build plot title
temp = strsplit(datapath,'/',/extract)
date = temp[4]
measurement = temp[5]
plot_title=strcompress(sample+'   '+date+'   '+measurement+'   channel '+string(channel))
plot_title='!I '+plot_title              ; define plot title and make it small

OutFileName = strcompress(sample+'   '+date+'   '+measurement+'   channel '+string(channel)+' analyzed (Excess QP).ps')
SummaryFileName = strcompress(sample+'   '+date+'   '+measurement+'   channel '+string(channel)+' summary (Excess QP).ps')

; set up plot
set_plot,'ps'
device,/color,encapsulated=0
loadct,4
!P.FONT = 0
device,/helv,/isolatin1      
device,font_size=12,/inches,xsize=7.5,ysize=9,xoffset=.5,yoffset=1
device,filename=outpath+OutFileName

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data = read_ascii(strcompress(datapath+'excess_TiN_QP_Fits_ch'+string(fix(channel))+'.dat',/remove_all))   ; read saved fits
data = data.field1

    ; fit: y(t) = (1-exp[-(t-t0)/Tau]) * [A*exp[-(t-t0)/Lambda] + D*A*[Lambda/(Lambda-Kappa)]*(exp[-(t-t0)/Lamda] - exp[-(t-t0)/Kappa]) + C*D*A*exp[-(t-t0)/Kappa]]

  ; p[0] = time offset of peak rise from 0 : t0 = p[0]-100
  ; p[1] = pulse rise time = Tau
  ; p[2] = peak height from thermal QP = A
  ; p[3] = fall time for membrane temp. = Lambda
  ; p[4] = ratio of energy in QP to energy in phonons = D
  ; p[5] = TiN QP lifetime = Kappa
  ; p[6] = amplification of the non-thermal QP signal with respect to the thermal signal = C

offset = (data[0,*])[*]   ; time offset of photon hit
rise = (data[1,*])[*]     ; pulse rise time = Tau
A = (data[2,*])[*]        ; peak height from thermal QP = A 
lambda = (data[3,*])[*]   ; fall time for membrane temp. = Lambda
D = (data[4,*])[*]        ; ratio of energy in QP to energy in phonons = D
kappa = (data[5,*])[*]    ; TiN QP lifetime = Kappa
C = (data[6,*])[*]        ; amplification of the non-thermal QP signal with respect to the thermal signal = C
chisq = (data[7,*])[*]    ; chi-squared     

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!p.multi=[0,1,2]    ; 2 diagrams per page

; histograms:
   ; d = histogram (x, min=a, max=b, binsize=c)  ==>> histogram of the array 'x': d[i] = number of elements in 'x' that have a value between c*i and c*(i+1)
   ;                                                  'a' = starting value of the first bin / 'b' = end value of the last bin ('a' and 'b' are not realy necessary)
   ; bs=1.0  ; bin size
   ; hist = histogram(fheight, MIN = 0, MAX = 220, BINSIZE = bs) ; y-values for the histogram plot
   ; bins = FINDGEN(N_ELEMENTS(hist))*bs                         ; x-values for the histogram plot  
   ; plot, bins, hist, psym=10, /ystyle, yr=[0,max(hist)*1.1], xtitle='Fast Fall Pulse Height (degrees)- all pulses, bin = 1', title=plot_title
   ;                      or (better locking & more simple):
   ; cghistoplot, x, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='xxx', title=plot_title, charsize=1
   
  ;  offset  rise  A  lambda  D  kappa  C  chisq

  filter = where( chisq LT MaxChisq and rise GT MinRiseTime)
 ; filter = where( chisq LT 35 and kappa GT 10 and lambda LT 800 and rise LT 500 and rise GT 0.2 and C LT 300)

  E = A[filter]+(D[filter]*A[filter])
 
 ; bs=1
 ; hist = E[0:100]
 ; cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='E - all membrane hits', title=plot_title, charsize=1;, xr=[0,500]
 

bs=1
hist = A+(D*A)
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='calculated photon energy - all membrane hits', title=plot_title, charsize=1;, xr=[0,100]

   bs= 2
   hist = A[filter]+(D[filter]*A[filter])
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='calculated photon energy - filtered further', title=plot_title, charsize=1;, xr=[0,200]
   

bs=1
hist = A
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='thermal pulse height - all membrane hits', title=plot_title, charsize=1;, xr=[0,0.5]

   bs=1
   hist = A[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='thermal pulse height - filtered further', title=plot_title, charsize=1;, xr=[3.8,4.2]
  

bs=0.05
hist = D
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='D: ratio QP to phonons  - all membrane hits', title=plot_title, charsize=1;, xr=[0,3]

   bs=0.01
   hist = D[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='D: ratio QP to phonons - filtered further', title=plot_title, charsize=1, xr=[0,2]
   
   
bs=1
hist = C
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='C: QP amp factor  - all membrane hits', title=plot_title, charsize=1;, xr=[0,50]

   bs=1
   hist = C[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='C: QP amp factor - filtered further', title=plot_title, charsize=1, xr=[0,100]
  
   
bs=1
hist = rise
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='rise time in uSec - all membrane hits', title=plot_title, charsize=1;, xr=[0,100]

   bs=1
   hist = rise[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='rise time in uSec - filtered further', title=plot_title, charsize=1;, xr=[0,100] 


bs=1.0
hist = lambda
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='thermal Lambda - all membrane hits', title=plot_title, charsize=1;, xr=[0,600]

   bs=1.0
   hist = lambda[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='thermal Lambda - filtered further', title=plot_title, charsize=1;, xr=[0,700]


bs=0.2
hist = kappa
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='QP liftetime Kappa - all membrane hits', title=plot_title, charsize=1;, xr=[0,100]

   bs=0.1
   hist = kappa[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='QP liftetime Kappa - filtered further', title=plot_title, charsize=1;, xr=[0,50]


bs=0.1
hist = chisq
cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='chi-squared - all membrane hits', title=plot_title, charsize=1;, xr=[0,100]

   bs=0.1
   hist = chisq[filter]
   cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='chi-squared - filtered further', title=plot_title, charsize=1;, xr=[0,30]

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

   ; psym values:
   ; 1 = +
   ; 2 = *
   ; 3 = .
   ; 4 = diamond
   ; 5 =  triangle
   ; 6 = square
   ; 7 = x
   ; 8 = user defined (usersym)
   ; 10 = histogram mode
   
if (scatter EQ 1) then begin    ;  offset  rise  A  lambda  D  kappa  C  chisq
!p.multi=[0,2,2]  ; 4 diagrams per page

plot, rise[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='Rise Time (!9m!3s)', ytitle='calculated photon energy', title=plot_title
plot, rise[filter], lambda[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(lambda[filter])*1.1], xtitle='Rise Time (!9m!3s)', ytitle='Thermal Fall Time (!9m!3s)', title=plot_title
plot, rise[filter], kappa[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(kappa[filter])*1.1], xtitle='Rise Time (!9m!3s)', ytitle='non-Thermal Fall Time (!9m!3s)', title=plot_title
plot, rise[filter], A[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(A[filter])*1.1], xtitle='Rise Time (!9m!3s)', ytitle='A', title=plot_title
plot, rise[filter], D[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(D[filter])*1.1], xtitle='Rise Time (!9m!3s)', ytitle='D', title=plot_title
plot, rise[filter], C[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,max(D[filter])*1.1], xtitle='Rise Time (!9m!3s)', ytitle='C', title=plot_title
plot, rise[filter], chisq[filter], psym=3, xr=[0,max(rise[filter])*1.1], yr=[0,5], xtitle='Rise Time (!9m!3s)', ytitle='Chi-squared', title=plot_title

plot, A[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(A[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='C', ytitle='calculated photon energy', title=plot_title
plot, A[filter], lambda[filter], psym=3, xr=[0,max(A[filter])*1.1], yr=[0,max(lambda[filter])*1.1], xtitle='Thermal Pulse Height', ytitle='Lambda (!9m!3s)', title=plot_title
plot, A[filter], kappa[filter], psym=3, xr=[0,max(A[filter])*1.1], yr=[0,max(kappa[filter])*1.1], xtitle='Thermal Pulse Height', ytitle='kappa (!9m!3s)', title=plot_title
plot, A[filter], D[filter], psym=3, xr=[0,max(A[filter])*1.1], yr=[0,max(D[filter])*1.1], xtitle='Thermal Pulse Height', ytitle='D', title=plot_title
plot, A[filter], C[filter], psym=3, xr=[0,max(A[filter])*1.1], yr=[0,max(C[filter])*1.1], xtitle='Thermal Pulse Height', ytitle='C', title=plot_title
plot, A[filter], chisq[filter], psym=3, xr=[0,max(A[filter])*1.1], yr=[0,5], xtitle='Thermal Pulse Height', ytitle='Chi-squared', title=plot_title

plot, lambda[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(lambda[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='Lambda (!9m!3s)', ytitle='calculated photon energy', title=plot_title
plot, lambda[filter], kappa[filter], psym=3, xr=[0,max(lambda[filter])*1.1], yr=[0,max(kappa[filter])*1.1], xtitle='Lambda (!9m!3s)', ytitle='kappa (!9m!3s)', title=plot_title
plot, lambda[filter], D[filter], psym=3, xr=[0,max(lambda[filter])*1.1], yr=[0,max(D[filter])*1.1], xtitle='Lambda (!9m!3s)', ytitle='D', title=plot_title
plot, lambda[filter], C[filter], psym=3, xr=[0,max(lambda[filter])*1.1], yr=[0,max(C[filter])*1.1], xtitle='Lambda (!9m!3s)', ytitle='C', title=plot_title
plot, lambda[filter], chisq[filter], psym=3, xr=[0,max(lambda[filter])*1.1], yr=[0,5], xtitle='Lambda (!9m!3s)', ytitle='Chi-squared', title=plot_title

plot, D[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(D[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='D', ytitle='calculated photon energy', title=plot_title
plot, D[filter], kappa[filter], psym=3, xr=[0,max(D[filter])*1.1], yr=[0,max(kappa[filter])*1.1], xtitle='D', ytitle='kappa (!9m!3s)', title=plot_title
plot, D[filter], C[filter], psym=3, xr=[0,max(D[filter])*1.1], yr=[0,max(C[filter])*1.1], xtitle='D', ytitle='C (!9m!3s)', title=plot_title
plot, D[filter], chisq[filter], psym=3, xr=[0,max(D[filter])*1.1], yr=[0,5], xtitle='D', ytitle='Chi-squared', title=plot_title

plot, kappa[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(kappa[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='kappa (!9m!3s)', ytitle='calculated photon energy', title=plot_title
plot, kappa[filter], C[filter], psym=3, xr=[0,max(kappa[filter])*1.1], yr=[0,max(C[filter])*1.1], xtitle='kappa (!9m!3s)', ytitle='C', title=plot_title
plot, kappa[filter], chisq[filter], psym=3, xr=[0,max(kappa[filter])*1.1], yr=[0,5], xtitle='kappa (!9m!3s)', ytitle='Chi-squared', title=plot_title

plot, C[filter], ((D[filter]*A[filter])+A[filter]), psym=3, xr=[0,max(C[filter])*1.1], yr=[0,max(((D[filter]*A[filter])+A[filter]))*1.1], xtitle='C (!9m!3s)', ytitle='calculated photon energy', title=plot_title
plot, C[filter], chisq[filter], psym=3, xr=[0,max(C[filter])*1.1], yr=[0,5], xtitle='C (!9m!3s)', ytitle='Chi-squared', title=plot_title

endif

device,/close

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;=====================================================================================================================================================================================================

if (summary EQ 1) then begin 

; set up plot
set_plot,'ps'
device,/color,encapsulated=0
loadct,4
!P.FONT = 0
device,/helv,/isolatin1      
device,font_size=12,/inches,xsize=7.5,ysize=9,xoffset=.5,yoffset=1
device,filename=outpath+SummaryFileName
temp = FINDGEN(17) * (!PI*2/16.) &  USERSYM, COS(temp), SIN(temp), /FILL ; defines filled circles as plot symbols
!p.multi=[0,1,2]    ; 2 diagrams per page

  
  ; >>>>>  plot photon energy determined from excess-QP-modell fit with fixed parameters 
  bs= 0.2
  hist = A[filter]+(D[filter]*A[filter])
  cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='calculated photon energy from membane hits', title=plot_title, charsize=1, xr=[0,1.1*max(hist)]
  aaa = strcompress('min. rise time = ' + string(MinRiseTime,format='(F6.1)') + ' !9m!3s')
  bbb = strcompress('max. Chi Sq = ' + string(MaxChisq,format='(F6.2)'))
  al_legend, /top, /left, [aaa,bbb], charsize=0.6 


  ; >>>>>  fit Gauss distribution to histogram to calculate R
  E = A[filter]+(D[filter]*A[filter])
  E1 = 0.0D  & E2 = 0.0D                            ; initiate E1 & E2 = distribution center peak 1 and peak 2
  sigma1 = 0.0D  & sigma2 = 0.0D                    ; sigmas for peak 1 and peak 2
     GaussFit = Gauss_fit_function(E, E1, sigma1, Amp1, E2, sigma2, Amp2, NBins)    
  ; scale from degrees to eV, used Fe55 numbers: 59.45% @ 5898.75 eV + 30.09% @ 5887.65 eV ( = 89.54% @ 5893.2) + 10.46% @ 6490.45 eV
  scale = 597.25 / ( E2-E1)                    ; scales values in degrees to values in eV
  shift = 5893.2 - E1 * 597.25 / (E2-E1)       ; shifts the rescaled eV values so that the first peak position E1 lies at 5893.2 eV  
     index = dindgen(max(E)*1.2*100)/100
     x = index * scale + shift                         
  y = histogram( E ,NBINS=Nbins ,LOCATIONS=BinValues)
  BinValues = BinValues + (Binvalues[1] - BinValues[0]) / 2.0        ; shift BinValues from 'start of bin' to 'center of bin'
  BinValues = BinValues * scale + shift           ; scale BinValues (in degrees) to eV 
  plot, BinValues, y, psym = 10, line = 0, thick = 2, color = 0, xr=[4500,7000], xtitle = 'photon energy [eV]', ytitle = 'no. of events', title = plot_title
  oplot, x, Gaussfit, color = 150, thick = 3
  Gauss1 = Amp1 * exp(-( (index-E1)^2 / (2*sigma1^2 ) ))   &   Gauss2 = Amp2 * exp(-( (index-E2)^2 / (2*sigma2^2 ) ))
  R1 = 5898.75 / (2.3548 * sigma1*scale)   &   R2 = 6490.45 / (2.3548 * sigma2*scale)
  R1leg = strcompress('R from 1. peak: ' + string(R1,format='(F6.1)'))
  R2leg = strcompress('R from 2. peak: ' + string(R2,format='(F6.1)'))
  al_legend, /top, /left, [R1leg, R2leg], charsize=0.6 


  ; >>>>>  plot maximum peak height of all membrane hits
  MembranePulseFile = strcompress(datapath+'ch'+string(fix(channel))+'_membrane_'+pulsename,/remove_all)  ; file name of membrane-hits data file
  N = long(((FILE_INFO(MembranePulseFile)).size)/(2*8.0*Npts))                                            ; every pulse is a set of 2 double (8 byte) variable (Pha,+ Amp) -->> number of pulses
  Pha = dblarr(Npts)
  Amp = dblarr(Npts)
  index = dindgen(Npts)  ; index for all measurement points
  index1 = dindgen(Npts/2-100)*1.25-125.0  ; index scaled and shifted: we measure with 800 kHz = 800.000 data point per Sec. = 1 data point every 1.25 uSec.
  openr, 1, MembranePulseFile  ; open membrane-hits data file to read
  maxPulse = dblarr(N)  &  maxPulse[*] = 0.0     ; pulse max in degrees
   for i=0L,N-1L do begin
    readu,1,Pha  ; read next already scaled phase pulse from file
    readu,1,Amp  ; read next already scaled amplitude pulse from file
    maxPulse[i] = max(Pha[Npts/2-100:Npts/2+100]*180.0/!Pi)    
    if (i EQ 3) then Pha_3 = Pha   ; save phase data of 3. membrane pulse
    if (i EQ 6) then Pha_6 = Pha   ; save phase data of 6. membrane pulse
    if (i EQ 9) then Pha_9 = Pha   ; save phase data of 9. membrane pulse
    if (i EQ 12) then Pha_12 = Pha   ; save phase data of 12. membrane pulse
    if (i EQ 15) then Pha_15 = Pha   ; save phase data of 15. membrane pulse
   endfor
  close,1
  bs = 1
  hist = maxPulse
  cghistoplot, hist, binsize=bs, datacolorname='red' ,POLYCOLOR='red', /fill, ytitle=' ', xtitle='just max. peak hight membrane hits - no fitting', title=plot_title, charsize=1, xr=[0,1.1*max(hist)]
  
  
  ; >>>>> plot 2 single Gauss fits to the histogram
  plot, BinValues, y, psym = 10, line = 0, thick = 2, color = 0, xr=[4500,7000], xtitle = 'photon energy [eV]', ytitle = 'no. of events', title = plot_title
  index = dindgen(max(E)*1.2*100)/100
  Gauss1 = Amp1 * exp(-( (index-E1)^2 / (2*sigma1^2 ) ))   &   Gauss2 = Amp2 * exp(-( (index-E2)^2 / (2*sigma2^2 ) ))
  oplot, x, Gauss1, color = 50, thick = 3
  oplot, x, Gauss2, color = 100, thick = 3
  E1Leg = strcompress('center 1. peak: ' + string( E1*((1+C[0]*D[0])/(1+D[0])),format='(F8.2)') + ' degrees, scaled to 5893.20 eV')    ; E1*((1+C[0]*D[0])/(1+D[0]) scales E1 = 'pulse max if C would be 0'  to  'max pulse if there would be no rise time'
  E2leg = strcompress('center 2. peak: ' + string( E2*((1+C[0]*D[0])/(1+D[0])),format='(F8.2)') + ' degrees, scaled to 6490.45 eV')
  Sig1leg = strcompress('sigma 1. peak: ' + string( sigma1*scale ,format='(F6.1)'))
  Sig2leg = strcompress('sigma 2. peak: ' + string( sigma2*scale ,format='(F6.1)'))
  ratio = strcompress('peak ratio : ' + string(Amp1/Amp2,format='(F4.1)') + ' : 1')
  saturation = strcompress('energy corresponding to 0 degrees phase shift : ' + string(shift,format='(F8.0)') + ' eV')
  al_legend, /top, /left, [E1leg, E2leg, Sig1leg, Sig2leg, ratio, saturation], charsize=0.6 
  index = dindgen(Npts)
  

  ; >>>>>  plot IQ-loop and amplitude-over-frequency for the frequency sweep with the same temp. and attenuation then the pulse file 
  iqsweepfit = strcompress(outpath + string(fix(temperature)) +'-'+ string(fix(StartRes)) +'-'+ string(fix(atten)) +'.ps',/remove_all)   ; file containing fit to frequqency sweep
  iqsweep = strcompress(datapath+string(fix(temperature)) +'-'+ string(fix(StartRes)) +'-'+ string(fix(atten)) +'.swp',/remove_all)      ; name of measured sweep data file
  lines = file_lines(iqsweep) - 7 ; number of points in frequency sweep, 7 lines of header are substracted
  ; We need the IQ-loop data and some data from the file header of the IQ-seeep-file.
  openr,1,iqsweep                    ; open IQ-loop file to read
  readf,1,fr1,fspan1,fsteps1,atten1  ; unused variables
  readf,1,fr2,fspan2,fsteps2,atten2  ; unused variables
  readf,1,ts,te     ; sweep start and end temperature
  readf,1,Iz1,Izsd1 ; origin of the IQ plain channel 1
  readf,1,Qz1,Qzsd1 ; origin of the IQ plain channel 1
  readf,1,Iz2,Izsd2 ; origin of the IQ plain channel 2
  readf,1,Qz2,Qzsd2 ; origin of the IQ plain channel 2
  loopdata = dblarr(5,lines) ; IQ-loop data
  readf,1,loopdata           ; read IQ-loop data
  close,1
  ; if the IQ-sweep was not measured as up-sweep but as down-sweep it's necessary to rearange loopdata:
  m = lines/2-1
  if( loopdata[0,0] GT loopdata[0,10] ) then begin ; if sweep was downsweep then ....
    print,'converting downsweep data to upsweep data'
    loopdata[0,0:m] = reverse(loopdata[0,0:m],2)
    loopdata[1,0:m] = reverse(loopdata[1,0:m],2)
    loopdata[2,0:m] = reverse(loopdata[2,0:m],2)
    loopdata[3,0:m] = reverse(loopdata[3,0:m],2)
    loopdata[4,0:m] = reverse(loopdata[4,0:m],2)  ; reverses the order in the array 'loopdata'
    loopdata[0,m+1:lines-1] = reverse(loopdata[0,m+1:lines-1],2)
    loopdata[1,m+1:lines-1] = reverse(loopdata[1,m+1:lines-1],2)
    loopdata[2,m+1:lines-1] = reverse(loopdata[2,m+1:lines-1],2)
    loopdata[3,m+1:lines-1] = reverse(loopdata[3,m+1:lines-1],2)
    loopdata[4,m+1:lines-1] = reverse(loopdata[4,m+1:lines-1],2)    
  endif
  if (channel EQ 1) then begin           ; channel 1 or channel 2
    loopdata_single = dblarr(5,lines/2)  &  loopdata_single[*,*] = loopdata[*,0:m]
    Izero = Iz1  &  Qzero = Qz1 
  endif else begin 
    loopdata_single = dblarr(5,lines/2)  &  loopdata_single[*,*] = loopdata[*,m+1:lines-1]
    Izero = Iz2  &  Qzero = Qz2  
  endelse
  I = (loopdata_single[1,*])[*]-Izero  ; extract I from data and substract origin shift
  Q = (loopdata_single[3,*])[*]-Qzero  ; extract Q from data and substract origin shift
  x = (loopdata_single[0,*])[*] * 1.0d9       ; measurement frequencies around rough f0
  y = dcomplex(I[*], Q[*]) ; complex measured IQ data
  iqfitname = datapath+'series-ps.dat' ; file name of the fit to the IQ sweep, file saved by the resfit program
  ; read data from fit to frequency sweep:
  loopfit = read_ascii(iqfitname) ; to read iqfitname like this will give a structure loopfit.field01, loopfit.field02, ...
  loopfit = loopfit.field01       ; get rid of unnecessary clutter in loopfit
    xc = loopfit[16,channel-1] + Izero       ; calculate coordinates of resonator loop center for both channels 
    yc = loopfit[17,channel-1] + Qzero
    f0 = loopfit[5,channel-1]           ; resonance frequency
    Q = loopfit[3,channel-1]            ; total Q
    Q_i = loopfit[9,channel-1]
    Q_c = loopfit[8,channel-1]
    aleak = loopfit[10,channel-1]    ;  amplitude of leakage
    ph1 = loopfit[11,channel-1]      ;  phase shift of leakage
    da = loopfit[12,channel-1]       ;  variation of carrier amplitude
    ang1 = loopfit[13,channel-1]     ;  rotation angle of data
    Igain = loopfit[14,channel-1]    ;  Gain of I channel
    Qgain = loopfit[15,channel-1]    ;  Gain of Q channel
    Ioff = loopfit[16,channel-1]     ;  Offset of I channel
    Qoff = loopfit[17,channel-1]     ;  Offset of Q channel
    db = loopfit[18,channel-1]       ; variation of carrier amplitude quadratic in frequency; not used at the moment  
  plot, double(y), imaginary(y), psym=8, title=strcompress(string(plot_title) + ' atten = ' + string(fix(atten)) +', T from '+string(ts*1000.0) + ' mK to ' + string(te*1000.0) + ' mK'), symsize=.2, /ynozero
  residx = where( min( (x-f0)^2 ) EQ (x-f0)^2 ) ; p[1]: f0 according to best fit / residx = where x-f0 minimal = measured frequency closest to f0
  plots, double(y[residx]), imaginary(y[residx]), psym=8, color=200 ; adds user-defined symbol at the f0-point in the IQ-loop-diagram
      ; calculate the fit (s21) to the measured IQ-loop using the function from 'RESDIFF':
        dx = (x - f0) / f0   ; shifts resonance to x=0
        s21a = (dcomplex(0,2.0*Q*dx)) / (dcomplex(1,0) + dcomplex(0,2.0*Q*dx))   ; resonator loop function
        s21a = s21a - dcomplex(.5,0)  ; shifts the loop so that its center lies at the IQ-origin
        s21b = dcomplex(da*dx + db*dx*dx,0) + s21a + aleak*dcomplex(1.0-cos(dx*ph1),-sin(dx*ph1)) ; adds the linear carrier amplitude variation and leakage
        ; scale, rotate and offset:
        Ix1 = double(s21b)*Igain     ; I-component of the function, scaled with the fitted gain of the I-channel
        Qx1 = imaginary(s21b)*Qgain  ; Q-component of the function, scaled with the fitted gain of the Q-channel
        nI1 = ((Ix1*cos(ang1) + Qx1*sin(ang1)))[*]  ; due to feedline resonances, the IQ-loops usually are rotated; this adds a rotation to the function
        nQ1 = ((-Ix1*sin(ang1) + Qx1*cos(ang1)))[*]
        nI1 = nI1 + Ioff  ; adds the fitted offset to I and Q
        nQ1 = nQ1 + Qoff
        s21 = dcomplex(nI1,nQ1) ; recombine I and Q to a complex value
  oplot, double(s21), imaginary(s21), color=150, thick=1 ; add s21 as line to the plot
  plots, Ioff, Qoff, psym=1, color=100  ; add a '+' to the plot at the position of the IQ-offset = center of IQ-loop
  ; plot measured and fitted resonance amplitude:
    amp = sqrt(double(y)^2 + imaginary(y)^2)   ; = abs(y) = measured resonance amplitude
    amp = amp/max(amp)  ; normalized
    amp = 20.0*alog10(amp)  ; in dB
    fitmag = sqrt(double(s21)^2 + imaginary(s21)^2)
    plot, x/1d9, amp, psym=8, /ynozero, /xstyle, symsize=.2, xtitle='Frequency (GHz)', ytitle='S!L21!N(dB)'  ; plot measured resonance amplitude (normalized and in dB) over f
    oplot, x/1d9, 20.0*alog10(fitmag/max(fitmag)), color=150                                                 ; add to plot: resonance amplitude (normalized and in dB) according to fit
    plots, x[residx]/1d9, amp[residx], psym=8, color=200                                                     ; add to plot: user defined symbol for measurement value at f0
    al_legend, /bottom, /right, [strcompress('Q = ' + string(Q, format='(F9.0)')), strcompress('Q!Lc!N = ' + string(Q_c, format='(F9.0)')), strcompress('Q!Li!N = ' + string(Q_i, format='(F9.0)')), strcompress('f!L0!N = ' + string(f0/1d9, format='(F9.6)') + ' GHz')], SPACING=1.5
  
  
  ; >>>>>  plot calculated photon energy of every membarne hit as time series (time axes not linear - in fact just pulse number) 
  y = A[filter]+(D[filter]*A[filter])
  plot, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='pulse no.', ytitle='calculated photon energy - membrane hits', symsize=0.25,  color = 50, ygridstyle = 1, yticklen = 1, thick = 0.1, /nodata
  oplot, y, psym=8, symsize=0.25, color = 150
  plots, [0,8600], [190,190], color = 50, thick = 3
  
  
  ; >>>>>  plot A+B for all pulses, A,B: amplitude of exponential functions from simple, all-parameters-free double-exponential fit 
  if (channel EQ 1) then data = read_ascii(datapath+'DoubleExpFits_ch1.dat') else data = read_ascii(datapath+'DoubleExpFits_ch2.dat')   ; read fitting parameters from saved double exp. fits
  data = data.field1    ; p[1] = pulse rise time / p[2] = A = amplitude first exponential / ; p[3] = Lambda / p[4] = B = amplitude 2. exponential
  rise_all = (data[1,*])[*]      ; rise time
  A_all = (data[2,*])[*]         ; amplitude first exponential
  lambda_all = (data[3,*])[*]    ; fall time 1. exponential
  B_all = (data[4,*])[*]         ; amplitude second exponential
  bs=1.0
  hist = A_all + B_all
  cghistoplot, hist, binsize=bs, datacolorname='darkgreen' ,POLYCOLOR='darkgreen', /fill, ytitle=' ', xtitle='A + B (double exp fit) - all pulses', title=plot_title, charsize=1, xr=[0,450];, yr=[0,N/20]


  ; >>>>>  plot rise time of all identified membrane pulses (rise time from excess-QP-modell fit)
  bs=0.5
  hist = rise[filter]
  cghistoplot, hist, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='rise time - membrane hits (fixed fits) [uSec]', title=plot_title, charsize=1, xr=[0,150];, yr=[0,N/20]
  
  ; >>>>>  plot rise time of all pulses (rise time from excess-QP-modell fit)
  bs=0.5
  hist = rise_all
  cghistoplot, hist, binsize=bs, datacolorname='blue' ,POLYCOLOR='blue', /fill, ytitle=' ', xtitle='rise time - all pulses (double exp. fits) [uSec]', title=plot_title, charsize=1, xr=[0,150], yr=[0,N]
  

  ; >>>>>  plot one of the 2 fall times from the simple, all-parameters-free double exponential fit for all pulses
  bs=1.0
  hist = lambda_all
  cghistoplot, hist, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='fall time 1. exponential (double exp fit) - all pulses [uSec]', title=plot_title, charsize=1, xr=[0,1500]


 

;  ; >>>>>  plot one membrane-hit pulse (arbitrary: pulse no. 3)
;  plot, index1, Pha_3[Npts/2-100:Npts-100]*180.0/!Pi, /xstyle, /ystyle, psym=8, xtitle='Time (microseconds)', ytitle='Phase Pulse Height (degrees)', symsize=0.25, title = plot_title, color=50   ; plots the measured data
;  j = 3  &  fit = (1-exp(-(x+100.0-offset[j])/rise[j])) * ((A[j]+D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/lambda[j]) + (C[j]*D[j]*A[j]-D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/kappa[j]))
;  oplot, index1, fit, line=0, color=150        ; adds the fitted function to the measured data   
;   rt = strcompress('rise time = ' + string(rise[j],format='(F6.1)') + ' !9m!3s')   
;   fh = strcompress('A = ' + string(A[j],format='(F6.1)') + ' dg')     
;   ff = strcompress('lambda = ' + string(lambda[j],format='(F6.1)') + ' !9m!3s')
;   sh = strcompress('D = ' + string(D[j],format='(F6.2)'))
;   sf = strcompress('kappa = ' + string(kappa[j],format='(F6.1)') + ' !9m!3s')       
;   tf = strcompress('C = ' + string(C[j],format='(F6.2)'))
;   csq = strcompress('Chi Sq = ' + string(chisq[j],format='(F6.2)'))        
;  al_legend, /top, /right, [rt,fh,ff,sh,sf,tf,csq], charsize=0.8 ; adds legend to the plot
  
  
  !p.multi=[0,2,3]
  x = double(index1) 

  ; >>>>>  plot 4 membrane hit pulses (arbitrary: pulses no. 6, 9, 12 and 15) 
  y = Pha_6[Npts/2-100:Npts-100]*180.0/!Pi
  plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height (degrees)', symsize=0.25, title = plot_title, color=50    ; plots the measured data
  j = 6  &  fit = (1-exp(-(x+100.0-offset[j])/rise[j])) * ((A[j]+D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/lambda[j]) + (C[j]*D[j]*A[j]-D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/kappa[j]))
  oplot, index1, fit, line=0, color=150        ; adds the fitted function to the measured data   
   rt = strcompress('rise time = ' + string(rise[j],format='(F6.1)') + ' !9m!3s')   
   fh = strcompress('A = ' + string(A[j],format='(F6.1)') + ' dg')     
   ff = strcompress('lambda = ' + string(lambda[j],format='(F6.1)') + ' !9m!3s')
   sh = strcompress('D = ' + string(D[j],format='(F6.2)'))
   sf = strcompress('kappa = ' + string(kappa[j],format='(F6.1)') + ' !9m!3s')       
   tf = strcompress('C = ' + string(C[j],format='(F6.2)'))
   csq = strcompress('Chi Sq = ' + string(chisq[j],format='(F6.2)'))        
  al_legend, /top, /right, [rt,fh,ff,sh,sf,tf,csq], charsize=0.8 ; adds legend to the plot

  y = Pha_9[Npts/2-100:Npts-100]*180.0/!Pi
  plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height (degrees)', symsize=0.25, title = plot_title, color=50    ; plots the measured data
  j = 9  &  fit = (1-exp(-(x+100.0-offset[j])/rise[j])) * ((A[j]+D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/lambda[j]) + (C[j]*D[j]*A[j]-D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/kappa[j]))
  oplot, index1, fit, line=0, color=150        ; adds the fitted function to the measured data   
   rt = strcompress('rise time = ' + string(rise[j],format='(F6.1)') + ' !9m!3s')   
   fh = strcompress('A = ' + string(A[j],format='(F6.1)') + ' dg')     
   ff = strcompress('lambda = ' + string(lambda[j],format='(F6.1)') + ' !9m!3s')
   sh = strcompress('D = ' + string(D[j],format='(F6.2)'))
   sf = strcompress('kappa = ' + string(kappa[j],format='(F6.1)') + ' !9m!3s')       
   tf = strcompress('C = ' + string(C[j],format='(F6.2)'))
   csq = strcompress('Chi Sq = ' + string(chisq[j],format='(F6.2)'))        
  al_legend, /top, /right, [rt,fh,ff,sh,sf,tf,csq], charsize=0.8 ; adds legend to the plot

  y = Pha_12[Npts/2-100:Npts-100]*180.0/!Pi
  plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height (degrees)', symsize=0.25, title = plot_title, color=50    ; plots the measured data
  j = 12  &  fit = (1-exp(-(x+100.0-offset[j])/rise[j])) * ((A[j]+D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/lambda[j]) + (C[j]*D[j]*A[j]-D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/kappa[j]))
  oplot, index1, fit, line=0, color=150        ; adds the fitted function to the measured data   
   rt = strcompress('rise time = ' + string(rise[j],format='(F6.1)') + ' !9m!3s')   
   fh = strcompress('A = ' + string(A[j],format='(F6.1)') + ' dg')     
   ff = strcompress('lambda = ' + string(lambda[j],format='(F6.1)') + ' !9m!3s')
   sh = strcompress('D = ' + string(D[j],format='(F6.2)'))
   sf = strcompress('kappa = ' + string(kappa[j],format='(F6.1)') + ' !9m!3s')       
   tf = strcompress('C = ' + string(C[j],format='(F6.2)'))
   csq = strcompress('Chi Sq = ' + string(chisq[j],format='(F6.2)'))        
  al_legend, /top, /right, [rt,fh,ff,sh,sf,tf,csq], charsize=0.8 ; adds legend to the plot

  y = Pha_15[Npts/2-100:Npts-100]*180.0/!Pi
  plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height (degrees)', symsize=0.25, title = plot_title, color=50    ; plots the measured data
  j = 15  &  fit = (1-exp(-(x+100.0-offset[j])/rise[j])) * ((A[j]+D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/lambda[j]) + (C[j]*D[j]*A[j]-D[j]*A[j]*(lambda[j]/(lambda[j]-kappa[j])))*exp(-(x+100.0-offset[j])/kappa[j]))
  oplot, index1, fit, line=0, color=150        ; adds the fitted function to the measured data   
   rt = strcompress('rise time = ' + string(rise[j],format='(F6.1)') + ' !9m!3s')   
   fh = strcompress('A = ' + string(A[j],format='(F6.1)') + ' dg')     
   ff = strcompress('lambda = ' + string(lambda[j],format='(F6.1)') + ' !9m!3s')
   sh = strcompress('D = ' + string(D[j],format='(F6.2)'))
   sf = strcompress('kappa = ' + string(kappa[j],format='(F6.1)') + ' !9m!3s')       
   tf = strcompress('C = ' + string(C[j],format='(F6.2)'))
   csq = strcompress('Chi Sq = ' + string(chisq[j],format='(F6.2)'))        
  al_legend, /top, /right, [rt,fh,ff,sh,sf,tf,csq], charsize=0.8 ; adds legend to the plot


  ; >>>>>  plot 2 (hopefully, not checked) non/membrane-hit pulses
  pulsefile = datapath+pulsename ; pulse data file name, contains pulses of 2 resonators (channel 1 and channel 2)
  openr,1,pulsefile  ; opens pulse data file
  header = dblarr(14)
  readu,1,header  ; reads the header of the pulse file -->> seems to be too big (112 byte instead of 56 byte), thus the first 12 pulses are read as 'header' and ignored -->> unimportant
  Ix1 = intarr(Npts)
  Qx1 = intarr(Npts)
  Ix2 = intarr(Npts)
  Qx2 = intarr(Npts)
   xc1 = loopfit[16,0] + Iz1       ; calculate coordinates of resonator loop center for both channels 
   yc1 = loopfit[17,0] + Qz1
   xc2 = loopfit[16,1] + Iz2
   yc2 = loopfit[17,1] + Qz2
  j = 0
  for i=0,20 do begin    ; 20 tries to find 2 substrate pulses
    readu,1,Ix1
    readu,1,Qx1
    readu,1,Ix2
    readu,1,Qx2 ; read next pulse from file
    Ix1d = (double(Ix1)/32767.0)*0.2
    Qx1d = (double(Qx1)/32767.0)*0.2
    Ix2d = (double(Ix2)/32767.0)*0.2
    Qx2d = (double(Qx2)/32767.0)*0.2  ; scale I and Q
     if( i EQ 0 ) then begin   ; only for the first pulse: 
        Ix1m = mean([Ix1d[0:Npts/2-100],Ix1d[Npts-100:-1]])
        Qx1m = mean([Qx1d[0:Npts/2-100],Qx1d[Npts-100:-1]])
        Ix2m = mean([Ix2d[0:Npts/2-100],Ix2d[Npts-100:-1]])
        Qx2m = mean([Qx2d[0:Npts/2-100],Qx2d[Npts-100:-1]])  ; calculate the mean values for I & Q where there is no pulse
        rad1 = sqrt( double(Qx1m-yc1)^2 + double(Ix1m-xc1)^2 ) 
        rad2 = sqrt( double(Qx2m-yc2)^2 + double(Ix2m-xc2)^2 ) ; calculate resonance depth = amplitude where ther is no signal
      endif    
    ; linear fit the noise to remove the base line, noise defined as "from start to 100 points before pulse" and "last 100 measurement points"
    r1 = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [Ix1d[0:Npts/2-100],Ix1d[Npts-100:-1]])
    r2 = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [Qx1d[0:Npts/2-100],Qx1d[Npts-100:-1]])
    r3 = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [Ix2d[0:Npts/2-100],Ix2d[Npts-100:-1]])
    r4 = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [Qx2d[0:Npts/2-100],Qx2d[Npts-100:-1]])
    ; substract linear fit AND add the mean noise value (this will be substracted later)
    Ix1d = Ix1d + Ix1m - (index*r1[1] + r1[0])
    Qx1d = Qx1d + Qx1m - (index*r2[1] + r2[0])
    Ix2d = Ix2d + Ix2m - (index*r3[1] + r3[0])
    Qx2d = Qx2d + Qx2m - (index*r4[1] + r4[0])
    ; transform I/Q to phase pulses:
    P1 = atan( double(Qx1d-yc1), double(Ix1d-xc1) )
    P2 = atan( double(Qx2d-yc2), double(Ix2d-xc2) )  ; simple angle calculation
    m1 = moment(P1[0:200])
    m2 = moment(P2[0:200]) ; m1[1], m2[1] = variance of the phase noise of the first 200 measured points
    P1 = fixang(P1,/radians)
    P2 = fixang(P2,/radians) ; removes discontinuities around +-Pi in angle values
   ; subtract phase baseline (again)
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P1[0:Npts/2-100],P1[Npts-100:-1]]) ; linear fit to noise in phase data
    P1 = P1 - index*temp[1] - temp[0]  ; substract linear noise component
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P2[0:Npts/2-100],P2[Npts-100:-1]])
    P2 = P2 - index*temp[1] - temp[0]
    ; reverse signs for A and P (looks better)
    P1 = -P1       
    P2 = -P2
    ; skip event if phase baseline sub is bad = if difference between first data point and 2/3 of (Npts/2-100) is bigger than 10 times the noise variance
    if( abs(P1[0]-P1[(Npts-200)/3]) GT 10.0*sqrt(m1[1]) ) then continue
    if( abs(P2[0]-P2[(Npts-200)/3]) GT 10.0*sqrt(m2[1]) ) then continue
    if(j EQ 2) then continue
  if (channel EQ 1) then begin
    if(max(P1[Npts/2-100:Npts/2+100]*180.0/!Pi) LT 40) then continue
    y = P1[Npts/2-100:Npts-100]*180.0/!Pi
    plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 1 (degrees)', symsize=0.25, color=200
    j++
   endif else begin
    if(max(P2[Npts/2-100:Npts/2+100]*180.0/!Pi) LT 40) then continue
    y = P2[Npts/2-100:Npts-100]*180.0/!Pi
    plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 2 (degrees)', symsize=0.25, color=200
    j++
   endelse
  endfor
  close, 1


device,/close
endif

print,' '
print,'  created plots in ',outpath+OutFileName
if (summary EQ 1) then print,' created summary file ',outpath+SummaryFileName
print,' '
print,'  done'
print,' '

end
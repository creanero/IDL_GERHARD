;*********************************************************************************************
;*  extract membrane hits from measurement.pro        Gerhard Ulbricht, January 2014         *
;*  based on  pulsefit.pro by Ben Mazin, June, 2013                                          * 
;*                                                                                           *
;*  This programm will read TKID pulse data from a pulse measurement file, add data from     *
;*  an IQ sweep file, bo a basic double exponential fit to all the pulses, identify          *
;*  membrane hits by their rise time and save these membrane pulses in new files             *
;*  "ch1/2_membrane_pulse_data.dat"                                                          *
;*                                                                                           *
;*                    <><><>      adjust before running:      <><><>                         *
;*                    <><><>    information in 'to_fit.txt'   <><><>                         *
;*                    <><><>  variable MembraneFilterTime1/2  <><><>                         *
;*                    <><><>    variable MinPulseHeight1/2    <><><>                         *
;*              <><><>  variable MinPulseTime1/2 or MinFallTime1/2  <><><>                   *
;*                  <><><>  desired membrane hit filter method  <><><>                       *
;*                    <><><>       variable ChisqCut/2        <><><>                         *
;*                    <><><>   all pulses or only 30.000?     <><><>                         *
;*                                                                                           *
;*********************************************************************************************

Function AtanWrap, P, index, Npts
; use to fix discontinuities in atan function intead of fixang
 
 P_return = P
 StDev = moment(P[0:(0.35*Npts)])  &  StDEv = sqrt(StDEv[1])    ; component 1 of procedure 'moment' = variance of a collection = normalized sum of squared deviations from mean value >> standard deviation = sqrt(variance) 
 PhaseMean = mean(P[0:(0.35*Npts)])
 
 ; decide: is a positive phase value (= phase rotation in 'wrong' direction) just noise or in fact a big rotation in right direction?
 temp = where((P[index] GT PhaseMean+2*StDev and P[index+1] GT PhaseMean+2*StDev and P[index+2] GT PhaseMean+2*StDev and P[index+3] GT PhaseMean+2*StDev) OR (P[index] GT PhaseMean+8*StDev), /Null) 
    ; all positive phase values that are in fact just big negative phase values >> they have to be shifted down by 360 degrees
 if temp NE !Null then P_return[temp] = P_return[temp]-2*!Pi
 return, P_return  

end

; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FUNCTION FittingFunction, p, X=x, Y=y, ERR=err

  ; used variables:
  ;  p: (input from outside) 6 fitting parameters p[0] - p[5]
  ;  x: (input from outside) times for each data point
  ;  y: (input from outside) measured data points
  ;  err: (input from outside) measurement error, ignored at the moment
  ;  f: function that should be fitted to data
  ;  dev: deviation between data point and fitted funktion divided by error

   
  ; general formula for an exponential decay with a final rise time:     (1-exp[-t/Tau]) * A*exp[-t/Lambda]
  ; Tau : pulse rise time (after Tau the pulse is at 63% of its theoretical max height, 86% after 2 Tau)
  ; 
  ; goal: fit every pulse with a generic, basic double exponential decay with a finite rise time
  ; doube exponential decay: f(t) = A*exp[-(t-t0)/Lambda] + B*exp[-(t-t0)/Kappa]
  ; t0 : where the pulse begins to rise
  ; Lambda, Kappa : decay times
  ; A, B : amplitudes
  ; 
  ; rise time added:
  ; f(t) = (1-exp[-(t-t0)/Tau]) * (A*exp[-(t-t0)/Lambda] + B*exp[-(t-t0)/Kappa])
  ; Tau: rise time


  ; p[0] = time offset of peak rise from 0 = t0
  ; p[1] = pulse rise time = Tau
  ; p[2] = amplitude first exponential = A
  ; p[3] = time constant 1. exponential = Lambda
  ; p[4] = amplitude 2. exponential = B
  ; p[5] = time constant 2. exponential = Kappa


  f = (1-exp(-(x+100.0-p[0])/p[1])) * (p[2]*exp(-(x+100.0-p[0])/p[3]) + p[4]*exp(-(x+100.0-p[0])/p[5]))

  k = where( f LT 0.0 ) ; set f(x) = 0 for x < 0 
  f[k]=0.0
  
  dev = y
  dev[*] = 0.0
  dev = (y-f)/err
  return, dev

end
; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

PRO FitIt, datapoints,index,p,fit,chisq, fitstatus
  ; calls the fit itself with the function defined in 'FittingFunction' and performs a complete least square fit on one pulse

  ; used variables:
  ;  datapoints = y: (input from outside) measured data
  ;  index = x: index, runs from 0 to Npts, used to number all measurement points
  ;  P: 6 fitting parameters p[0] - p[5]
  ;  fit: y-values for the best fit
  ;  chisq: Chi-squared calculated from the best fit
  ;  fitstatus: status variable of the fit
  ;  err: relative measurement error, ignored (set to 1.0) at the moment
  ;  parinfo: array (6 x structure with 6 elements) necessary for the fitting algorithm, contains starting values and borders for all 6 fitting parameters
  ;  fa1: structure containing x,y and err; input for fitting algorithm
  ;  bestnorm: summed squared and weighted deviations between best fit and data points
  ;  covar: covariance matrix for best fit parameters  
  ;  perror: formal 1-sigma errors in each parameter, inaccurate as long as all error are set to 1.0  
  ;  dof: degrees of freedom = Npts - 6    

  x = double(index)  ; make sure x and y are double precision
  y = double(datapoints)
  err = replicate(1.0, (size(x))[1] )  ; set all errors to 1.0, effectively ignoring them

    parinfo = replicate({value:0.D, fixed:0, step:0, limited:[1,1], limits:[0.D,0.D],mpmaxstep:0.D}, 6)
    ; ...value: start value, ...fixed: to fit or not to fit, ...step: step size for numerical drivation (0 for automatic)
    ; ...limited: bounded to lower/upper bound or not, ...limits: upper/lower bounds, ...mpmaxstep: max parameter change per iteration (0: no max)
    parinfo[0].value = 100.0
    parinfo[0].limits  = [80.0,120.0]  ; peak time offset from 0 uSec between -20 and +20 (fitting function substracts 100 uSec from p[0])
    
    parinfo[1].value = 1.0
    parinfo[1].limits  = [0.01,200.0] ; rise time between 0.01 and 60 uSec

    parinfo[2].value = 75.0
    parinfo[2].limits  = [1.0,360.0] ; A between 0.01 and 260 degrees

    parinfo[3].value = 200.0
    parinfo[3].limits  = [1,3000.0] ; Lambda between 5 and 3000 uSec
  
    parinfo[4].value = 75.0
    parinfo[4].limits  = [1.0,360.0] ; B between 0.01 and 100
  
    parinfo[5].value = 200.0
    parinfo[5].limits  = [1.0,3000.0] ; Kappa between 20 and 250 uSec
    

  ; p[0] = time offset of peak rise from 0 = t0
  ; p[1] = pulse rise time = Tau
  ; p[2] = amplitude first exponential = A
  ; p[3] = time constant 1. exponential = Lambda
  ; p[4] = amplitude 2. exponential = B
  ; p[5] = time constant 2. exponential = Kappa


  ; do the optimization:
  fa1 = {X:x, Y:y, ERR:err}
  bestnorm=0.0
  covar=0.0
  perror=0.0 ; first set all outputs of the mpfit procedure to 0

  p = mpfit('FittingFunction', FUNCTARGS=fa1, BESTNORM=bestnorm, COVAR=covar, DOF=dof, PARINFO=parinfo, PERROR=perror, AUTODERIVATIVE=1, FTOL=1D-18, XTOL=1D-18, GTOL=1D-18, FASTNORM=0, STATUS=status, /QUIET) 
  ; does a complete least square fit of the function 'FittingFunction' to the data points in fa1, returns p[0] - p[5] as best fit parameters
  fitstatus = status

   ; print,'Status of Fit= ',status
    ;print,bestnorm
    ;DOF     = N_ELEMENTS(X) - N_ELEMENTS(PARMS) ; deg of freedom
    ;PCERROR = PERROR * SQRT(BESTNORM / DOF)   ; scaled uncertainties
    ;print,dof,bestnorm   ; lines for debugging
        ; fit = p[4]*(1.0-exp(-(x-100.0+p[2])/p[0]))*exp(-(x-100.0+p[2])/p[5]) + p[3]*(1.0-exp(-(x-100.0+p[2])/p[0]))*exp(-(x-100.0+p[2])/p[1])
  fit = (1-exp(-(x+100.0-p[0])/p[1])) * (p[2]*exp(-(x+100.0-p[0])/p[3]) + p[4]*exp(-(x+100.0-p[0])/p[5]))
  k = where( fit LT 0.0 )
  fit[k]=0.0
  chisq = bestnorm/dof

end

; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FUNCTION SingleExponential, p_exp, X=x, Y=y, ERR=err

  ; used variables:
  ;  p_exp: (input from outside) 4 fitting parameters p_exp[0] - p_exp[3]
  ;  x: (input from outside) times for each data point
  ;  y: (input from outside) measured data points
  ;  err: (input from outside) measurement error, ignored at the moment
  ;  f_exp: function that should be fitted to data
  ;  dev: deviation between data point and fitted funktion divided by error

  ; simple single exponential decay with final rise time:    (1-exp[-t-t0/Tau]) * A*exp[-t-t0/Lambda]
  ; used to test how many pulses are better (or as well) fitted with a single then a double exponential decay
  ; Tau : pulse rise time
  ; t0 : where the pulse begins to rise
  ; Lambda: exsp. decay time
  ; A: exp. amplitude

  ; p_exp[0] = time offset of peak rise from 0 = t0
  ; p_exp[1] = pulse rise time = Tau
  ; p_exp[2] = amplitude = A
  ; p_exp[3] = time constant = Lambda

  f_exp = (1-exp(-(x+100.0-p_exp[0])/p_exp[1])) * p_exp[2]*exp(-(x+100.0-p_exp[0])/p_exp[3])
    
  k = where( f_exp LT 0.0 ) ; set f_exp(x) = 0 for x < 0 
  f_exp[k]=0.0
  
  dev = y
  dev[*] = 0.0
  dev = (y-f_exp)/err
  return, dev

end
; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

PRO FitSingleExp, datapoints,index,p_exp,fit_exp,chisq_exp, fitstatus_exp
  ; calls a single exponential fit to see how many pulses are better described with a single fit then with a double exponential

  ; used variables:
  ;  datapoints = y: (input from outside) measured data
  ;  index = x: index, runs from 0 to Npts, used to number all measurement points
  ;  P_exp: 6 fitting parameters p[0] - p[5]
  ;  fit_exp: y-values for the best fit
  ;  chisq_exp: Chi-squared calculated from the best fit
  ;  fitstatus_exp: status variable of the fit
  ;  err: relative measurement error, ignored (set to 1.0) at the moment
  ;  parinfo: array (6 x structure with 6 elements) necessary for the fitting algorithm, contains starting values and borders for all 6 fitting parameters
  ;  fa1: structure containing x,y and err; input for fitting algorithm
  ;  bestnorm: summed squared and weighted deviations between best fit and data points
  ;  covar: covariance matrix for best fit parameters  
  ;  perror: formal 1-sigma errors in each parameter, inaccurate as long as all error are set to 1.0  
  ;  dof: degrees of freedom = Npts - 6    

  x = double(index)  ; make sure x and y are double precision
  y = double(datapoints)
  err = replicate(1.0, (size(x))[1] )  ; set all errors to 1.0, effectively ignoring them

    parinfo = replicate({value:0.D, fixed:0, step:0, limited:[1,1], limits:[0.D,0.D],mpmaxstep:0.D}, 4)
    ; ...value: start value, ...fixed: to fit or not to fit, ...step: step size for numerical drivation (0 for automatic)
    ; ...limited: bounded to lower/upper bound or not, ...limits: upper/lower bounds, ...mpmaxstep: max parameter change per iteration (0: no max)
    parinfo[0].value = 100.0
    parinfo[0].limits  = [80.0,120.0]  ; peak time offset from 0 uSec between -20 and +20 (fitting function substracts 100 uSec from p[0])
    
    parinfo[1].value = 1.0
    parinfo[1].limits  = [0.01,60.0] ; rise time between 0.01 and 60 uSec

    parinfo[2].value = 75.0
    parinfo[2].limits  = [10.0,360.0] ; A between 0.01 and 260 degrees

    parinfo[3].value = 450.0
    parinfo[3].limits  = [10,3000.0] ; Lambda between 5 and 3000 uSec
  
  ; p_exp[0] = time offset of peak rise from 0 = t0
  ; p_exp[1] = pulse rise time = Tau
  ; p_exp[2] = amplitude = A
  ; p_exp[3] = time constant = Lambda


  ; do the optimization:
  fa1 = {X:x, Y:y, ERR:err}
  bestnorm=0.0
  covar=0.0
  perror=0.0 ; first set all outputs of the mpfit procedure to 0

  p_exp = mpfit('SingleExponential', FUNCTARGS=fa1, BESTNORM=bestnorm, COVAR=covar, DOF=dof, PARINFO=parinfo, PERROR=perror, AUTODERIVATIVE=1, FTOL=1D-18, XTOL=1D-18, GTOL=1D-18, FASTNORM=0, STATUS=status, /QUIET) 
  ; does a complete least square fit of the function 'FittingFunction' to the data points in fa1, returns p[0] - p[5] as best fit parameters
  fitstatus_exp = status

  fit_exp = (1-exp(-(x+100.0-p_exp[0])/p_exp[1])) * p_exp[2]*exp(-(x+100.0-p_exp[0])/p_exp[3])
  k = where( fit_exp LT 0.0 )
  fit_exp[k]=0.0
  chisq_exp = bestnorm/dof

end

; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

PRO FitPulses, datapath, outpath, pulsename, StartT, atten, res, Npts, sample
  ; procedure to call the least square fit for every single pulse

  ; used variables:
  ;  datapath: path where to find the data file
  ;  outpath: path where to save output file
  ;  pulsename: file containing pulse data to be analyzed
  ;  pulsefile: pulse data file name + path
  ;  StartT: measurement temperature in mK
  ;  atten: attenuator setting for pulse file in dB
  ;  res: resonator group as numbered by the IQ sweep routine
  ;  Npts: number of points in each pulse
  ;  sampe: sample name
  ;  N: number of pulses
  ;  iqsweep: file name for the correct IQ sweep data for the pulses
  ;  iqfitname: file name of the fit to the IQ sweep
  ;  Iz1, Iz2, Qz1, Qz2: origin of the IQ plain (channels 1 & 2), possibly shifted due to DC offsets
  ;  loopfit: array which contains values of the fit to the IQ-sweep loop
  ;  xc1, xc2, yc1, yc2: coordinates of both resonator IQ-loop centers
  ;  temp: temporary variable, only used local
  ;  header: header of the pulse data file
  ;  max1, max2: maximum peak heights in channel 1 / 2 in degrees
  ;  Ix1, Qx1, Ix2, Qx2: measured I and Q values for channel 1/2 for every data point 
  ;  Ix1d, Qx1d, Ix2d, Qx2d: Ix1, Qx1, Ix2, Qx2 in double precision and scaled
  ;  index: index, runs from 0 to Npts, used to number all measurement points
  ;  index1: index scaled and shifted for a correct representation of the measurement times
  ;  Ix1m, Qx1m, Ix2m, Qx2m: mean value for I/Q-noise where ther is no pulse
  ;  rad1, rad2: radius of the IQ-loop of channel 1/2 = resonance depth where there is no pulse
  ;  r1 - r4: linear fits to the noise in Ix1, Qx1, Ix2, Qx2 = base lines
  ;  P1, P2: measured resonance phase, transormed from measured I/Q values
  ;  A1, A2: measured resonance amplitude, transormed from measured I/Q values
  ;  m1[1], m2[1]: variance of the phase noise of the first 200 measured points
  ;  P: 6 pulse fitting parameters p[0] - p[5]
  ;  fit: y-values for the best fit
  ;  chisq: Chi-squared calculated from the best fit
  ;  fitstatus: status variable of the fit
  ;  p_exp, fit_exp, chisq_exp, fitstatus_exp: p, fit, chisq and fitstatus for the single exponential fit
  ;  rt, fh, ff, sh, sf, csq: strings only used for plot legend
  ;  scip1, scip2: scip channel 1/2 pulse because of bad baseline fit
  ;  baselinedecline1, baselinedecline2: counter for declined pulses due to bad baseline fit
  ;  pulseheightdecline1, pulseheightdecline2: counter for pulses declined because of too small height
  ;  MembraneFilterTime1, 2: time in uSec to determine if a pulse was a membrane hit or not (membrane hit if risetime > MembraneFilterTime)
  ;  MinPulseHeight1, 2: ignore all pulses that are smaller then this in phase-degrees
  ;  MembranePulseFile1, 2: file to save the measurement data of identified membrane pulses into  
  ;  MembranePulseHeader1, 2: file to save the parameters used to identify channel 1, 2 membrane pulses
  ;  SingleChisq1, 2: number of membrane pulses where the single exponential fit gives a chi-squared of less then ChisqCut1, 2
  ;  DoubleChisq1, 2: number of membrane pulses where the double exponential fit gives a chi-squared of less then ChisqCut1, 2
  ;  ChisqCut1, 2: fit is considered a good fit if Chi-squared below that value
  

  MembraneFilterTime1 = 3.0D    ; smallest rise time (in uSec) in order to identify a pulse as membrane hit, channel 1
  MembraneFilterTime2 = 3.0D  ; smallest rise time (in uSec) in order to identify a pulse as membrane hit, channel 2
  MinPulseHeight1 = 40.0      ; channel 1: only fit if pulse is greater than this   (AND only identify as membrane hit if pulse bigger then this MinPulseTime1 in uSec after photon hit)
  MinPulseHeight2 = 40.0      ; channel 2: only fit if pulse is greater than this   (AND only identify as membrane hit if pulse bigger then this MinPulseTime2 in uSec after photon hit)
  ; MinPulseTime1 = 40          ; 1 data point every 1.25 uSec, so 50 uSec are 50/1.25 = 40 points  
  ; MinPulseTime2 = 40          ; 1 data point every 1.25 uSec, so 50 uSec are 50/1.25 = 40 points  
  ; MinFallTime1 = 100          ; min. second fall time required for membrane hit channel 1
  ; MinFallTime2 = 100          ; min. second fall time required for membrane hit channel 2
  MinChisq1 = 200.0            ; min. Chisquared for membrane hit Ch1
  MinChisq2 = 200.0            ; min. Chisquared for membrane hit Ch1
  ChisqCut1 = 2.0             ; channel 1 : fit is considered a good fit if Chi-squared below that value (only used for counting, not filtering)
  ChisqCut2 = 2.0             ; channel 2 : fit is considered a good fit if Chi-squared below that value (only used for counting, not filtering)
  
  
  iqsweep = strcompress(datapath+string(fix(StartT)) +'-'+ string(fix(res)) +'-'+ string(fix(atten)) +'.swp',/remove_all) ; IQ-sweep file name
  pulsefile = datapath+pulsename ; pulse data file name, contains pulses of 2 resonators (channel 1 and channel 2)
  iqfitname = datapath+'series-ps.dat' ; file name of the fit to the IQ sweep, file saved by the resfit program
  MembranePulseFile1 = datapath+'ch1_membrane_'+pulsename   ; file to save the measurement data of channel 1 identified membrane pulses into 
  MembranePulseFile2 = datapath+'ch2_membrane_'+pulsename   ; file to save the measurement data of channel 2 identified membrane pulses into 
  MembranePulseHeader1 = datapath+'parameters_for_ch1_membrane_'+pulsename   ; file to save the parameters used to identify channel 1 membrane pulses
  MembranePulseHeader2 = datapath+'parameters_for_ch2_membrane_'+pulsename   ; file to save the parameters used to identify channel 2 membrane pulses
  
  fiterr1 = 0L  &  fiterr2 = 0L    ; set several counters to zero
  fitmaxiter1 = 0L  &  fitmaxiter2 = 0L
  fitsuccess1 = 0L  &  fitsuccess2 = 0L
  baselinedecline1 = 0L  &  baselinedecline2 = 0L
  pulseheightdecline1 = 0L  &  pulseheightdecline2 = 0L
  NumMembraneHits1 = 0L  &  NumMembraneHits2 = 0L
  SingleChisq1 = 0L  &  SingleChisq2 = 0L
  DoubleChisq1 = 0L  &  DoubleChisq2 = 0L

  ; We need the origin of the IQ plane in order to calculate phase and amplitude from I and Q. The origin is the measured value for no excitation
  ; and not exact 0 due to DC offsets. The layout of the IQ-seeep-file is: 7 lines of header, followed by the IQ sweep data in a table with 5 columns
  ; the header consists of: channel 1 sweep start frequqncy - sweep end frequency - frequqncy step size - attentuator setting
  ;                         channel 2 sweep start frequqncy - sweep end frequency - frequqncy step size - attentuator setting
  ;                         start temperature - end temperature
  ;                         channel 1: I-value for 'no exitation' - standard deviation (= noise) of that value
  ;                         channel 1: Q-value for 'no exitation' - standard deviation (= noise) of that value
  ;                         channel 2: I-value for 'no exitation' - standard deviation (= noise) of that value
  ;                         channel 2: Q-value for 'no exitation' - standard deviation (= noise) of that value
  ; to read the desired values Iz1, Iz2, Qz1 and Qz2 it's best to read the whole header and just ignore / not use all the other variables:
  openr,1,iqsweep
  readf,1,fstart1,fend1,fsteps1,atten1
  readf,1,fstart2,fend2,fsteps2,atten2
  readf,1,Tstart,Tend
  readf,1,Iz1,Izsd1
  readf,1,Qz1,Qzsd1
  readf,1,Iz2,Izsd2
  readf,1,Qz2,Qzsd2
  close,1

  ; we also need the center of the resonator loop according to its fit, this value can be found in the iqfitname file:
  loopfit = read_ascii(iqfitname) ; to read iqfitname like this will give a structure loopfit.field01, loopfit.field02, ...
  loopfit = loopfit.field01       ; get rid of unnecessary clutter in loopfit
  xc1 = loopfit[16,0] + Iz1       ; calculate coordinates of resonator loop center for both channels 
  yc1 = loopfit[17,0] + Qz1
  xc2 = loopfit[16,1] + Iz2
  yc2 = loopfit[17,1] + Qz2

  ; set up postscript output
  set_plot,'ps'
  device,/color,encapsulated=0
  loadct,4
  !p.multi=[0,2,3]
  !P.FONT = 0
  device,/helv,/isolatin1      
  device,font_size=12,/inches,xsize=7.5,ysize=9,xoffset=.5,yoffset=1
  device,filename=datapath+'BasicDoubleExpFits.ps'
  temp = FINDGEN(17) * (!PI*2/16.) &  USERSYM, COS(temp), SIN(temp), /FILL ; defines filled circles as plot symbols

  openr,1,pulsefile  ; opens pulse data file
  header = dblarr(14)
  readu,1,header  ; reads the header of the pulse file -->> seems to be too big (112 byte instead of 56 byte), thus the first 12 pulses are read as 'header' and ignored -->> unimportant

  N = long(((FILE_INFO(pulsefile)).size-56.0)/(4*2.0*Npts)) - 1  
  ; N = 30000   ; don't work with the complete file, only use the first 30.000 pulses instead
  ; the file has a 56 bit header & every measurement consist of 4 values (Ix1, Qx1, Ix2, Qx2) that are all integer = 2 byte long
  ; N = long(1000) ; just analayze the first xxx pulses
  print,'Reading',N,' pulses'

  ; initiate some variables used later
  max1 = dblarr(N)
  max2 = dblarr(N)
  max1[*] = 0.0
  max2[*] = 0.0
  E1 = dblarr(N)
  E2 = dblarr(N)
  E1[*] = -10.0
  E2[*] = -10.0    ; I later plot only possitive values of E, thus ignoring all E values that have not been changed because pulse was too small
  count1 = 0  &  count2 = 0  &  count3 = 0  &  count4 = 0  
  count5 = 0  &  count6 = 0  &  count7 = 0  &  count8 = 0
  count9 = 0  &  count10 = 0  & count11 = 0  & count12 = 0
  Ix1 = intarr(Npts)
  Qx1 = intarr(Npts)
  Ix2 = intarr(Npts)
  Qx2 = intarr(Npts)

  readu,1,Ix1
  readu,1,Qx1
  readu,1,Ix2
  readu,1,Qx2 ; read 1. pulse from file
  Ix1d = (double(Ix1)/32767.0)*0.2
  Qx1d = (double(Qx1)/32767.0)*0.2
  Ix2d = (double(Ix2)/32767.0)*0.2
  Qx2d = (double(Qx2)/32767.0)*0.2  ; scale Ix1-Qx2 and make them double precision
  ; the measurement instruments measure I and Q as 16-bit number (so, between -32767 and +32767) with a range of -0.2 V to +0.2 V

  index = dindgen(Npts)                    ; index for all measurement points
  index1 = dindgen(Npts/2-100)*1.25-125.0  ; index scaled and shifted
  ; we measure with 800 kHz = 800.000 data point per Sec. = 1 data point every 1.25 uSec.

  openw,2,datapath+'DoubleExpFits_ch1.dat'  ; generates and opens path+DoubleExpFits_ch1.dat to save fit data of channel 1 in -->> check later if MembraneFilterTime was set to a easonable value
  openw,3,datapath+'DoubleExpFits_ch2.dat'  ; generates and opens path+DoubleExpFits_ch2.dat to save fit data of channel 2 in -->> check later if MembraneFilterTime was set to a easonable value
  openw,4,MembranePulseFile1                ; generates and opens MembranePulseFile1 to save the measurement data of identified channel 1 membrane pulses into  
  openw,5,MembranePulseFile2                ; generates and opens MembranePulseFile2 to save the measurement data of identified channel 2 membrane pulses into  
  
  openw,6,MembranePulseHeader1              ; file to save the parameters used to identify channel 1 membrane pulses
  openw,7,MembranePulseHeader2              ; file to save the parameters used to identify channel 2 membrane pulses
  temp = 'parameters used in data file  ' + MembranePulseFile1 + '  to identify membrane hits on sample' + sample      ; save parameters for both MembranePulseFile containing used to identify memebrane hits
  printf, 6, temp
  temp = 'parameters used in data file  ' + MembranePulseFile2 + '  to identify membrane hits on sample' + sample      ; save parameters for both MembranePulseFile containing used to identify memebrane hits
  printf, 7, temp  
  temp = 'file format: 1. line = calculated phase signal, 2. line = calculated amplitude signal, both after baseline substraction'
  printf, 6, temp  & printf, 7, temp
  temp = 'original pulse data file:  ' + pulsefile  &  printf, 6, temp  & printf, 7, temp
  temp = 'channel 1'  &  printf, 6, temp
  temp = 'channel 2'  &  printf, 7, temp
  temp = 'used filter criteria: working base line substraction and fitting procedure'  &    printf, 6, temp  & printf, 7, temp
  temp = 'minimum pulse height in degrees:' + string(MinPulseHeight1)  &  printf, 6, temp
  temp = 'minimum pulse height in degrees:' + string(MinPulseHeight2)  &  printf, 7, temp
  temp = 'minimum pulse rise time in uSec:' + string(MembraneFilterTime1)  &  printf, 6, temp
  temp = 'minimum pulse rise time in uSec:' + string(MembraneFilterTime2)  &  printf, 7, temp
  temp =  MembranePulseFile1 + ' has no header and every single pulse consits of' + string(fix(Npts)) + '  points'  &    printf, 6, temp 
  temp =  MembranePulseFile2 + ' has no header and every single pulse consits of' + string(fix(Npts)) + '  points'  &    printf, 7, temp 
  temp =  'data was written and has to be read out as double-precision floating point'  &    printf, 6, temp   &   printf, 7, temp 
  
; ........................................................................................................................................
  
  ; start a loop with 1 iteration per pulse:
  for i=0L,N-1L do begin
    readu,1,Ix1
    readu,1,Qx1
    readu,1,Ix2
    readu,1,Qx2 ; read next pulse from file
    Ix1d = (double(Ix1)/32767.0)*0.2
    Qx1d = (double(Qx1)/32767.0)*0.2
    Ix2d = (double(Ix2)/32767.0)*0.2
    Qx2d = (double(Qx2)/32767.0)*0.2  ; scale I and Q
    
    ; debugging: plot Ix1d & Ix2d as read from the file
    ; plot, index, Ix1d[0:-1], /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Ix1d', symsize=0.25
    ; plot, index, Ix2d[0:-1], /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Ix2d', symsize=0.25
    
      if( i mod 50 EQ 0 ) then print,i,' /',N,'  membrane hits so far:',NumMembraneHits1,' /',NumMembraneHits2  ; print N for every 50th pulse and the number of identified membrane hits

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
    ; substract linear fit AND add the mean noise value to avoid shifting the loops around in the IQ plane by substratcing the linear fit 
    Ix1d = Ix1d + Ix1m - (index*r1[1] + r1[0])
    Qx1d = Qx1d + Qx1m - (index*r2[1] + r2[0])
    Ix2d = Ix2d + Ix2m - (index*r3[1] + r3[0])
    Qx2d = Qx2d + Qx2m - (index*r4[1] + r4[0])
    
    ; debugging: plot Ix1d after linear fit
    ; plot, index, Ix1d[0:-1], /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Ix1d', symsize=0.25

    ; transform I/Q to phase pulses:
    P1 = atan( double(Qx1d-yc1), double(Ix1d-xc1) )
    P2 = atan( double(Qx2d-yc2), double(Ix2d-xc2) )  ; simple angle calculation
    m1 = moment(P1[0:200])
    m2 = moment(P2[0:200]) ; m1[1], m2[1] = variance of the phase noise of the first 200 measured points
    P1_bevore_fixang = P1  &  P2_bevore_fixang = P2
    P1 = fixang(P1,/radians)
    P2 = fixang(P2,/radians) ; removes discontinuities around +-Pi in angle values

    ; transform I/Q to amplitude pulses (normalized to first pulse):
    A1 = sqrt( double(Qx1d-yc1)^2 + double(Ix1d-xc1)^2 ) / rad1
    A2 = sqrt( double(Qx2d-yc2)^2 + double(Ix2d-xc2)^2 ) / rad2
    
    ; debugging: plot phase signal before linear fit
    ; plot, index, P1[0:-1], /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Phase Pulse Height after calculation', symsize=0.25

    ; subtract phase baseline (again)
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P1[0:Npts/2-100],P1[Npts-100:-1]]) ; linear fit to noise in phase data
    P1 = P1 - index*temp[1] - temp[0]  ; substract linear noise component
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P2[0:Npts/2-100],P2[Npts-100:-1]])
    P2 = P2 - index*temp[1] - temp[0]
    ; subtract amplitude baseline (again)
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [A1[0:Npts/2-100],A1[Npts-100:-1]])
    A1 = A1 - index*temp[1] - temp[0]
    temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [A2[0:Npts/2-100],A2[Npts-100:-1]])
    A2 = A2 - index*temp[1] - temp[0]
    ; reverse signs for A and P (looks better)
    P1 = -P1       ; sometimes it's necesaary to comment these out (reason unclear)
    P2 = -P2
    A1 = -A1  &  A2 = -A2
    
    ; debugging: plot phase signal after linear fit
    ; plot, index, P1[0:-1]*180.0/!Pi, /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Phase Pulse Height ch1 (degrees)', symsize=0.25
    ; plot, index, P2[0:-1]*180.0/!Pi, /xstyle, /ystyle, psym=8, xtitle='Time (unscaled)', ytitle='Phase Pulse Height ch2 (degrees)', symsize=0.25
    
    ; skip event if phase baseline sub is bad = if difference between first data point and 2/3 of (Npts/2-100) is bigger than 10 times the noise variance
    scip1=0.0 & scip2=0.0
    
    if( abs(P1[0]-P1[(Npts-200)/3]) GT 10.0*sqrt(m1[1]) ) then begin
      print,'pulse',i,' in channel 1 has base line issues, trying AtanWrap'       
      P1 = AtanWrap(P1_bevore_fixang, index, Npts)       
        temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P1[0:Npts/2-100],P1[Npts-100:-1]]) ; linear fit to noise in phase data
        P1 = P1 - index*temp[1] - temp[0]                                                            ; substract linear noise component
      P1 = -P1      
      if( abs(P1[0]-P1[(Npts-200)/3]) GT 10.0*sqrt(m1[1]) ) then begin
        print,'pulse',i,' in channel 1 still has base line issues,  AtanWrap did not help'  &  print,''
        baselinedecline1++  
        scip1=1.0
        if(count9 LT 4) then begin
          y = P1[Npts/2-100:Npts-100]*180.0/!Pi
          plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase channel 1 (degrees), bad baseline fit', symsize=0.25, color=0               
          al_legend, /top, /left, ['bad base line'], charsize=0.8 
          count9++
        endif
      endif else begin 
        print,'AtanWrap got rid of base line issues'  &  print,''
        if(count10 LT 4) then begin
          y = P1[Npts/2-100:Npts-100]*180.0/!Pi
          plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase channel 1 (degrees), AtanWrap worked', symsize=0.25, color=0               
          al_legend, /top, /left, ['AtanWrap worked'], charsize=0.8 
          count10++
        endif
      endelse   
    endif

    if( abs(P2[0]-P2[(Npts-200)/3]) GT 10.0*sqrt(m2[1]) ) then begin
      print,'  pulse',i,' in channel 2 has base line issues, trying AtanWrap'       
      P2 = AtanWrap(P2_bevore_fixang, index, Npts)       
        temp = linfit( [index[0:Npts/2-100],index[Npts-100:-1]], [P2[0:Npts/2-100],P2[Npts-100:-1]]) ; linear fit to noise in phase data
        P2 = P2 - index*temp[1] - temp[0]                                                            ; substract linear noise component
      P2 = -P2      
      if( abs(P2[0]-P2[(Npts-200)/3]) GT 10.0*sqrt(m2[1]) ) then begin
        print,'  pulse',i,' in channel 2 still has base line issues,  AtanWrap did not help'  &  print,''
        baselinedecline2++  
        scip2=1.0
        if(count11 LT 4) then begin
          y = P2[Npts/2-100:Npts-100]*180.0/!Pi
          plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase channel 2 (degrees), bad baseline fit', symsize=0.25, color=0               
          al_legend, /top, /left, ['bad base line'], charsize=0.8 
          count11++
        endif
      endif else begin 
        print,'  AtanWrap got rid of base line issues'  &  print,''
        if(count12 LT 4) then begin
          y = P2[Npts/2-100:Npts-100]*180.0/!Pi
          plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase channel 2 (degrees), AtanWrap worked', symsize=0.25, color=0               
          al_legend, /top, /left, ['AtanWrap worked'], charsize=0.8 
          count12++
        endif
      endelse   
    endif

    ; get quick pulse max for each pulse i, only looking 100 points around peak trigger = +-125 uSec
    max1[i] = max(P1[Npts/2-100:Npts/2+100]*180.0/!Pi)
    max2[i] = max(P2[Npts/2-100:Npts/2+100]*180.0/!Pi)

  ; p[0] = time offset of peak rise from 0 = t0
  ; p[1] = pulse rise time = Tau
  ; p[2] = amplitude first exponential = A
  ; p[3] = time constant 1. exponential = Lambda
  ; p[4] = amplitude 2. exponential = B
  ; p[5] = time constant 2. exponential = Kappa   

    ; -- fit and save data for channel 1 --
    ; *************************************
    if (max1[i] GT MinPulseHeight1 AND scip1 EQ 0.0) then begin ; only fit a pulse if it is bigger then 20 degrees
      p=0
      fit=0
      chisq=0
      ; start procedure 'FitIt' to fit measurement data P1 (phase signal channel 1) between measurement point 'Npts/2-100' (100 points before pulse
      ; trigger) and 'Npts-100'; the corresponding x-values are index1 = equidistant steps between 0 and Npts-100 scaled to 'pulse trigger at 0' and 
      ; 'one point each 1.25 uSec'. Results will be returned in p (7 fitting parameters), fit (y values best fit) and chisq :
      FitIt, P1[Npts/2-100:Npts-100]*180.0/!Pi, index1, p, fit, chisq, fitstatus
      ; if one tries to fit to the last measurement point the fit can run into trouble as soon as the pulse starts a little late because then the data and
      ; the fit function would have different length.
      E1[i] = p[2] + p[4]  ; rough number for photon energy
      if (fitstatus LE 0) then fiterr1++     ; count the unsuccessfull fits
      if (fitstatus EQ 5) then fitmaxiter1++
      
      if (fitstatus GT 0 AND fitstatus NE 5) then begin    ; fit was a success -->> save the pulse data if it was a membrane hit
        fitsuccess1++
        
       ; if (p[1] GT MembraneFilterTime1 and P1[Npts/2+(p[0]-100)/1.25 + MinPulseTime1]*180.0/!Pi GT MinPulseHeight1) then begin   
                            ; identify membrane hit: (p[1] = pulse rise time) membrane hit if rise time bigger then MembraneFilterTime and Phase after MinPulseTime1 uSec still above MinPulseHeight
                            ; time of photon hit = Npts/2+t0 >> t0 = 100-p[0] in uSec = (100-p[0])/1.25 in P1-index (1 data point every 1.25 uSec)   (50/1.25=40)
                            
       ; if (p[1] GT MembraneFilterTime1 and p[5] GT MinFallTime1) then begin  
        if (p[1] GT MembraneFilterTime1 and chisq LT MinChisq1) then begin  
                            ; identify membrane hit: (p[1] = pulse rise time) membrane hit if rise time bigger then MembraneFilterTime and second exp. fall time above MinFallTime  >> alternative filter method
                            
          writeu, 4, P1     ; save channel 1 phase data of identified membrane pulse into the file 'MembranePulseFile'
          writeu, 4, A1     ; save channel 1 amplitude data of identified membrane pulse into the file 'MembranePulseFile'
          NumMembraneHits1++

          FitSingleExp, P1[Npts/2-100:Npts-100]*180.0/!Pi, index1, p_exp, fit_exp, chisq_exp, fitstatus_exp
                       ; for debugging: plot some single exp. pulses:
                       if( count6 LT 40 ) then begin
                        y = P1[Npts/2-100:Npts-100]*180.0/!Pi
                        plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 1 (degrees), membrane hit', symsize=0.25, color=50
                        oplot, index1, fit_exp, line=0, color=0
                        rt = strcompress('rise time = ' + string(p_exp[1],format='(F6.1)') + ' !9m!3s')
                        fh = strcompress('A = ' + string(p_exp[2],format='(F6.1)') + ' dg')
                        ff = strcompress('lambda = ' + string(p_exp[3],format='(F6.1)') + ' !9m!3s')
                        csq = strcompress('Chi Sq = ' + string(chisq_exp,format='(F6.2)'))
                        al_legend, /top, /right, [rt,fh,ff,csq], charsize=0.8 
                        count6++
                       endif
          if (chisq_exp LT ChisqCut1) then SingleChisq1++  ; count the number of good single exponential fits, membrane hits only
          if (chisq LT ChisqCut1 ) then DoubleChisq1++     ; count the number of good double exponential fits, membrane hits only

          if( count1 LT 40 ) then begin      ; plot the first 40 membrane pulses to BasicDoubleExpFits.ps
            y = P1[Npts/2-100:Npts-100]*180.0/!Pi
            plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 1 (degrees), membrane hit', symsize=0.25, color = 100       ; plots the measured data
            oplot, index1, fit, line=0, color=0  ; adds the fitted function to the measured data
            rt = strcompress('rise time = ' + string(p[1],format='(F6.1)') + ' !9m!3s')   
            fh = strcompress('A = ' + string(p[2],format='(F6.1)') + ' dg')     
            ff = strcompress('lambda = ' + string(p[3],format='(F6.1)') + ' !9m!3s')
            sh = strcompress('B = ' + string(p[4],format='(F6.3)'))
            sf = strcompress('kappa = ' + string(p[5],format='(F6.1)') + ' !9m!3s')       
            csq = strcompress('Chi Sq = ' + string(chisq,format='(F6.2)'))        
            al_legend, /top, /right, [rt,fh,ff,sh,sf,csq], charsize=0.8 ; adds legend to the plot
            count1++ 
          endif

          ; select an exemplary pulse and save it to file 'example_pulse.dat'
          ;if( count5 EQ 0.0 and chisq LT 1.0 ) then begin
          if( count5 EQ 0.0 ) then begin
           openw, 8, datapath+'example_pulse_double_exp_ch1.dat'
           printf, 8, P1, format='(7F9.3)'
           close, 8
           count5 = 1
           y = P1[Npts/2-100:Npts-100]*180.0/!Pi
           plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='example pulse ch1, membrane hit', symsize=0.25, color=0      ; plots the measured data
           oplot, index1, fit, line=0, color=150  ; adds the fitted function to the measured data
           print, 'example pulse ch1 saved'
          endif 
          
        endif else begin
          if( count2 LT 25 ) then begin      ; plot the first 20 not-membrane pulses to BasicDoubleExpFits.ps
            y = P1[Npts/2-100:Npts-100]*180.0/!Pi
            plot, index1, P1[Npts/2-100:Npts-100]*180.0/!Pi, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 1 (degrees), not-a-membrane hit', symsize=0.25, color = 80     ; plots the measured data
            oplot, index1, fit, line=0, color=0  ; adds the fitted function to the measured data
            rt = strcompress('rise time = ' + string(p[1],format='(F6.1)') + ' !9m!3s')   
            fh = strcompress('A = ' + string(p[2],format='(F6.1)') + ' dg')     
            ff = strcompress('lambda = ' + string(p[3],format='(F6.1)') + ' !9m!3s')
            sh = strcompress('B = ' + string(p[4],format='(F6.1)'))
            sf = strcompress('kappa = ' + string(p[5],format='(F6.1)') + ' !9m!3s')       
            csq = strcompress('Chi Sq = ' + string(chisq,format='(F6.2)'))        
            al_legend, /top, /right, [rt,fh,ff,sh,sf,csq], charsize=0.8 ; adds legend to the plot
            count2++
          endif  
          
        endelse  
      endif
      
      ; save to file path+'DoubleExpFits_ch1.dat': all 6 fitting parameters and chisq (for each pulse)
      printf, 2, p[0], p[1], p[2], p[3], p[4], p[5], chisq, format='(7F9.3)'  
      
    endif else pulseheightdecline1++
           
          
    ; -- fit and save data for channel 2 --
    ; *************************************   
    if (max2[i] GT MinPulseHeight2 AND scip2 EQ 0.0 ) then begin ; only fit a pulse if it is bigger then 20 degrees
      p=0
      fit=0
      chisq=0
      ; start procedure 'FitIt' to fit measurement data P2 (phase signal channel 2) between measurement point 'Npts/2-100' (100 points before pulse
      ; trigger) and 'Npts-100'; the corresponding x-values are index1 = equidistant steps between 0 and Npts-100 scaled to 'pulse trigger at 0' and 
      ; 'one point each 1.25 uSec'. Results will be returned in p (7 fitting parameters), fit (y values best fit) and chisq :
      FitIt, P2[Npts/2-100:Npts-100]*180.0/!Pi, index1, p, fit, chisq, fitstatus
      E2[i] = p[2] + p[4]  ; rough number for photon energy
      if (fitstatus LE 0) then fiterr2++     ; count the unsuccessfull fits
      if (fitstatus EQ 5) then fitmaxiter2++
  
      if (fitstatus GT 0 AND fitstatus NE 5) then begin    ; fit was a success -->> save the pulse data if it was a membrane hit
       fitsuccess2++
       
      ; if (p[1] GT MembraneFilterTime2 and P2[Npts/2+(p[0]-100)/1.25 + MinPulseTime2]*180.0/!Pi GT MinPulseHeight2) then begin   
                           ; identify membrane hit: (p[1] = pulse rise time) membrane hit if rise time bigger then MembraneFilterTime and Phase after MinPulseTime2 uSec still above MinPulseHeight
                           ; time of photon hit = Npts/2+t0 >> t0 = 100-p[0] in uSec = (100-p[0])/1.25 in P1-index (1 data point every 1.25 uSec)   (50/1.25=40)
                           
       ; if (p[1] GT MembraneFilterTime2 and p[5] GT MinFallTime2) then begin 
       if (p[1] GT MembraneFilterTime2 and chisq LT MinChisq2) then begin   
                            ; identify membrane hit: (p[1] = pulse rise time) membrane hit if rise time bigger then MembraneFilterTime and second exp. fall time above MinFallTime  >> alternative filter method
                           
         writeu, 5, P2     ; save channel 2 phase data of identified membrane pulse into the file 'MembranePulseFile'
         writeu, 5, A2     ; save channel 2 amplitude data of identified membrane pulse into the file 'MembranePulseFile'
         NumMembraneHits2++
         
         FitSingleExp, P2[Npts/2-100:Npts-100]*180.0/!Pi, index1, p_exp, fit_exp, chisq_exp, fitstatus_exp
                       ; for debugging: plot some single exp. pulses:
                       if( count7 LT 40 ) then begin     
                        y = P2[Npts/2-100:Npts-100]*180.0/!Pi
                        plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 2 (degrees), membrane hit', symsize=0.25, color=150
                        oplot, index1, fit_exp, line=0, color=0
                        rt = strcompress('rise time = ' + string(p_exp[1],format='(F6.1)') + ' !9m!3s')
                        fh = strcompress('A = ' + string(p_exp[2],format='(F6.1)') + ' dg')
                        ff = strcompress('lambda = ' + string(p_exp[3],format='(F6.1)') + ' !9m!3s')
                        csq = strcompress('Chi Sq = ' + string(chisq_exp,format='(F6.2)'))
                        al_legend, /top, /right, [rt,fh,ff,csq], charsize=0.8 
                        count7++
                       endif
          if (chisq_exp LT ChisqCut2) then SingleChisq2++  ; count the number of good single exponential fits, membrane hits only
          if (chisq LT ChisqCut2) then DoubleChisq2++     ; count the number of good double exponential fits, membrane hits only
         
         if( count3 LT 40 ) then begin      ; plot the first 40 membrane pulses to BasicDoubleExpFits.ps
           y = P2[Npts/2-100:Npts-100]*180.0/!Pi
           plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 2 (degrees), membrane hit', symsize=0.25, color = 200       ; plots the measured data
           oplot, index1, fit, line=0, color=0  ; adds the fitted function to the measured data
           rt = strcompress('rise time = ' + string(p[1],format='(F6.1)') + ' !9m!3s')   
           fh = strcompress('A = ' + string(p[2],format='(F6.1)') + ' dg')     
           ff = strcompress('lambda = ' + string(p[3],format='(F6.1)') + ' !9m!3s')
           sh = strcompress('B = ' + string(p[4],format='(F6.3)'))
           sf = strcompress('kappa = ' + string(p[5],format='(F6.1)') + ' !9m!3s')       
           csq = strcompress('Chi Sq = ' + string(chisq,format='(F6.2)'))        
           al_legend, /top, /right, [rt,fh,ff,sh,sf,csq], charsize=0.8 ; adds legend to the plot
           count3++ 
         endif
         
         ; select an exemplary pulse and save it to file
         if( count8 EQ 0.0 ) then begin
         ; if( count8 EQ 0.0 and p[1] GT 3.0 and chisq LT 1.5) then begin  
          openw, 9, datapath+'example_pulse_double_exp_ch2.dat'
          printf, 9, P2, format='(7F9.3)'
          close, 9
          count8 = 1
          y = P2[Npts/2-100:Npts-100]*180.0/!Pi
          plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='example pulse ch2, membrane hit', symsize=0.25, color=0      ; plots the measured data
          oplot, index1, fit, line=0, color=150  ; adds the fitted function to the measured data
          print, 'example pulse ch2 saved'
         endif             
         
       endif else begin
         if( count4 LT 25 ) then begin      ; plot the first 20 not-membrane pulses to BasicDoubleExpFits.ps
           y = P2[Npts/2-100:Npts-100]*180.0/!Pi
           plot, index1, y, /xstyle, /ystyle, psym=8, yr=[min(y),1.1*max(y)], xtitle='Time (microseconds)', ytitle='Phase Pulse Height channel 2 (degrees), not-a-membrane hit', symsize=0.25,  color = 125   ; plots the measured data
           oplot, index1, fit, line=0, color=0  ; adds the fitted function to the measured data
           rt = strcompress('rise time = ' + string(p[1],format='(F6.1)') + ' !9m!3s')   
           fh = strcompress('A = ' + string(p[2],format='(F6.1)') + ' dg')     
           ff = strcompress('lambda = ' + string(p[3],format='(F6.1)') + ' !9m!3s')
           sh = strcompress('B = ' + string(p[4],format='(F6.3)'))
           sf = strcompress('kappa = ' + string(p[5],format='(F6.1)') + ' !9m!3s')       
           csq = strcompress('Chi Sq = ' + string(chisq,format='(F6.2)'))        
           al_legend, /top, /right, [rt,fh,ff,sh,sf,csq], charsize=0.8 ; adds legend to the plot
           count4++
         endif  
         
       endelse  
     endif
      
      ; save to file path+'DoubleExpFits_ch2.dat': all 6 fitting parameters and chisq (for each pulse)
      printf, 3, p[0], p[1], p[2], p[3], p[4], p[5], chisq, format='(7F9.3)'
      
    endif else pulseheightdecline2++
 
     
  ; p[0] = time offset of peak rise from 0 = t0
  ; p[1] = pulse rise time = Tau
  ; p[2] = amplitude first exponential = A
  ; p[3] = time constant 1. exponential = Lambda
  ; p[4] = amplitude 2. exponential = B
  ; p[5] = time constant 2. exponential = Kappa   
  
  endfor
  
  ; plot max pulse heights to check for bunching
  bs=1.0
  cghistoplot, max1, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='max pulse height ch1', charsize=1, yrange=[0,N/20]
  cghistoplot, max2, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='max pulse height ch2', charsize=1, yrange=[0,N/20]
  
  ; plot rough energy value E1/2 to check for bunching
  bs=1.0
  cghistoplot, E1, binsize=bs, datacolorname='navy' ,POLYCOLOR='navy', /fill, ytitle=' ', xtitle='rough photon energy (A+B) ch1', charsize=1, xr=[0,400]
  cghistoplot, E2, binsize=bs, datacolorname='orange' ,POLYCOLOR='orange', /fill, ytitle=' ', xtitle='rough photon energy (A+B) ch2', charsize=1, xr=[0,400]
     
  
; ........................................................................................................................................

  close,1  &  close,2  &  close,3 
  close,4  &  close, 5             ; close all open files for saved pulses and fitting parameters
  device,/close ; close PostScript file (it contains the plots)
  
  ; save the number of successfull fits, .... in the MembranePulseHeader-File
  temp = ' successfull pulse fits:'+string(fitsuccess1)+'  --  max iterations reached:'+string(fitmaxiter1)+'  --  unsuccessfull fits:'+string(fiterr1)  &  printf, 6, temp 
  temp = ' successfull pulse fits:'+string(fitsuccess2)+'  --  max iterations reached:'+string(fitmaxiter2)+'  --  unsuccessfull fits:'+string(fiterr2)  &  printf, 7, temp  
  temp = ' pulses in channel 1 declined due to bad phase baseline substraction:'+string(baselinedecline1)   &  printf, 6, temp 
  temp = ' pulses in channel 2 declined due to bad phase baseline substraction:'+string(baselinedecline2)   &  printf, 7, temp   
  temp = ' pulses in channel 1 not fitted as they are too small:'+string(pulseheightdecline1)   &  printf, 6, temp 
  temp = ' pulses in channel 2 not fitted as they are too small:'+string(pulseheightdecline2)   &  printf, 7, temp 
  temp = ' identified and saved membrane hits:'+string(NumMembraneHits1)   &  printf, 6, temp   
  temp = ' identified and saved membrane hits:'+string(NumMembraneHits2)   &  printf, 7, temp   
  temp = ' single exponential fits for channel 1 with Chi-squared <'+string(fix(ChisqCut1))+' : '+string(SingleChisq1)  &  printf, 6, temp  
  temp = ' single exponential fits for channel 2 with Chi-squared <'+string(fix(ChisqCut2))+' : '+string(SingleChisq2)  &  printf, 7, temp    
  temp = ' double exponential fits for channel 1 with Chi-squared <'+string(fix(ChisqCut1))+' : '+string(DoubleChisq1)  &  printf, 6, temp    
  temp = ' double exponential fits for channel 2 with Chi-squared <'+string(fix(ChisqCut2))+' : '+string(DoubleChisq2)  &  printf, 7, temp   
  close, 6  &  close,7
  
  print,' '
  print,' channel 1: successfull pulse fits:',fitsuccess1,'  --  max iterations reached:',fitmaxiter1,'  --  unsuccessfull fits:',fiterr1  ; usual values: 8000 pulses, 1300 successfull fits, 0-2 max iterations reached, 0 unseuccessfull
  print,' channel 2: successfull pulse fits:',fitsuccess2,'  --  max iterations reached:',fitmaxiter2,'  --  unsuccessfull fits:',fiterr2
  print,' '
  print,' pulses in channel 1 declined due to bad phase baseline substraction:',baselinedecline1  ; usual value: 10 (broken analog readout: 150)
  print,' pulses in channel 2 declined due to bad phase baseline substraction:',baselinedecline2
  print,' '
  print,' pulses in channel 1 not fitted as they are below',MinPulseHeight1,' degrees:',pulseheightdecline1  
  print,' pulses in channel 2 not fitted as they are below',MinPulseHeight2,' degrees:',pulseheightdecline2
  print,' '
  print,' identified and saved membrane hits in channel 1:',NumMembraneHits1
  print,' identified and saved membrane hits in channel 2:',NumMembraneHits2
  print,' '
  print,' single exponential fits for channel 1 with Chi-squared < ',fix(ChisqCut1),' : ',SingleChisq1  
  print,' double exponential fits for channel 1 with Chi-squared < ',fix(ChisqCut1),' : ',DoubleChisq1  
  print,' '
  print,' single exponential fits for channel 2 with Chi-squared < ',fix(ChisqCut2),' : ',SingleChisq2    
  print,' double exponential fits for channel 2 with Chi-squared < ',fix(ChisqCut2),' : ',DoubleChisq2    
  print,' '  
  print,' created ',datapath+'DoubleExpFits_ch1.dat : --  channel 1 fitting parameters'
  print,' created ',datapath+'DoubleExpFits_ch2.dat : --  channel 2 fitting parameters'
  Print,''
  print,' created ',MembranePulseFile1 +' : --  measurement data of identified channel 1 membrane pulses'
  print,' created ',MembranePulseFile2 +' : --  measurement data of identified channel 2 membrane pulses'
  print,'' 
  print,' created ',MembranePulseHeader1 +' : --  parameters used to identify channel 1 membrane pulses'
  print,' created ',MembranePulseHeader2 +' : --  parameters used to identify channel 2 membrane pulses'
  print,''
  print,' created ',datapath+'example_pulse_double_exp_ch1.dat : --  example pulse channel 1'
  print,' created ',datapath+'example_pulse_double_exp_ch2.dat : --  example pulse channel 2'  
  print,''
  print,' created ',datapath+'BasicDoubleExpFits.ps : -- plots of the first membrane and not-membrane pulses'
  print,' '

; p.s.: whoever finds any spelling errors in the comments may keep them ;)
end

; ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

; main program

; used variables:
;  path: path to the master file
;  masterfile: name of the file that defines which data to fit
;  datapath:  path where to find the data file
;  outpath: path where to save output file
;  pulsename: file containing pulse data to be analyzed
;  temp: temporary variable, only used local
;  StartRes: resonator group as numbered by the IQ sweep routine
;  StartAtten: attenuator setting for pulse file in dB
;  StartT: measurement temperature in mK
;  bg: do a background subtraction - temp for high temp sweep (?), 0 for none, -1 for phase fit only (xxxx)
;  loss: attenuation before device, DilutionFridge = 45, ADR = 35  (XXXX)
;  Npts: number of points in each pulse
;  iqsweep: file containing the data of the frequency sweep = resonator IQ loop
;  loopdata: measurements of the resonators IQ loop
;  iqsweepfit: PS-file containing the fit to  iqsweep
;  ts, te: IQ-sweep start and end temperature
;  Iz1, Iz2, Qz1, Qz2: origin of the IQ plain (channels 1 & 2), possibly shifted due to DC offsets
;  lines: the numbers of measurement points in the IQ sweep
;  m: lines /2 -1
;  data1a, Iz1a, Qz1a, Iz2a, Qz2a: only used for a different background substraction of the IQ sweep data, unused at the moment
;  sample: sample name



; define a masterfile that defines which data to work with
path = '/home/gerhard/' ; where to find the file 'to_fit.txt'
cd,path
masterfile = 'to_fit.txt'
line=string(200)

; read all necessary information from masterfile
openr,1,strcompress(path+masterfile)
readf,1,line & datapath = (strsplit(line,/EXTRACT))[0] ; path where to find the data file
readf,1,line & outpath = (strsplit(line,/EXTRACT))[0]  ; path where to save output file
readf,1,line & pulsename = (strsplit(line,/EXTRACT))[0] ; file with pulse data to be analyzed
readf,1,line & temp = (strsplit(line,'resnum=',/EXTRACT))[0] & reads,temp,StartRes ; resonator group as numbered by the IQ sweep routine
readf,1,line & temp = (strsplit(line,'atten=',/EXTRACT))[0] & reads,temp,StartAtten ; attenuator setting for pulse file in dB
readf,1,line & temp = (strsplit(line,'T=',/EXTRACT))[0] & reads,temp,StartT ; measurement temperature in mK
readf,1,line & temp = (strsplit(line,'sub=',/EXTRACT))[0] & reads,temp,bg ; do a background subtraction or not
readf,1,line & temp = (strsplit(line,'loss=',/EXTRACT))[0] & reads,temp,loss ; attenuation before device, DilutionFridge = 45, ADR = 35
readf,1,line & temp = (strsplit(line,'Np=',/EXTRACT))[0] & reads,temp,Npts  ; number of points in each pulse
readf,1,line & sample = (strsplit(line,/EXTRACT))[0]   ; read sample name
close,1

; fit sweep data = resonator IQ loop, but only if the fit file doesn't already exist  
iqsweepfit = strcompress(outpath + string(fix(StartT)) +'-'+ string(fix(StartRes)) +'-'+ string(fix(StartAtten)) +'.ps',/remove_all)
if ((FILE_INFO(iqsweepfit)).EXISTS EQ 0) then begin
  iqsweep = strcompress(datapath+string(fix(StartT)) +'-'+ string(fix(StartRes)) +'-'+ string(fix(StartAtten)) +'.swp',/remove_all) ; name of sweep data file

  ; check if IQ-sweep-file is present, if not abort 
  lines = file_lines(iqsweep) - 7 ; 7 lines of header are substracted
  if( lines LE 100 ) then begin
    print,'File' + iqsweep + 'is missing or corrupt.'
    stop 
  endif

  ; We need the IQ-loop data and some data from the file header of the IQ-seeep-file. The header has 7 lines: 
  ; 		channel 1 sweep start frequqncy - sweep end frequency - frequency step size - attentuator setting
  ;     channel 2 sweep start frequqncy - sweep end frequency - frequqncy step size - attentuator setting
  ;     start temperature - end temperature
  ;     channel 1: I-value for 'no exitation' - standard deviation (= noise) of that value
  ;     channel 1: Q-value for 'no exitation' - standard deviation (= noise) of that value
  ;     channel 2: I-value for 'no exitation' - standard deviation (= noise) of that value
  ;     channel 2: Q-value for 'no exitation' - standard deviation (= noise) of that value
  ; followed by the IQ sweep data in a table with 5 columns:
  ; loopdata [0]: frequency
  ; loopdata [1]: measured I
  ; loopdata [2]: error in I
  ; loopdata [3]: measured Q
  ; loopdata [4]: error in Q
  openr,1,iqsweep  ; open IQ-loop file to read
  readf,1,fr1,fspan1,fsteps1,atten1  ; unused variables
  readf,1,fr2,fspan2,fsteps2,atten2  ; unused variables
  readf,1,ts,te  ; sweep start and end temperature
  readf,1,Iz1,Izsd1 ; origin of the IQ plain channel 1
  readf,1,Qz1,Qzsd1 ; origin of the IQ plain channel 1
  readf,1,Iz2,Izsd2 ; origin of the IQ plain channel 2
  readf,1,Qz2,Qzsd2 ; origin of the IQ plain channel 2
  loopdata = dblarr(5,lines) ; IQ-loop data
  readf,1,loopdata ; read IQ-loop data
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

  data1a = dblarr(5,lines)
  data1a[*,*] = 0.0  ; unused at the moment

  ; set error values in the IQ loop data to effectively zero:
  loopdata[2, * ] = 1d-6  &  loopdata[4, * ] = 1d-6

  ; postscript output setup
  set_plot,'ps'
  device,/color,encapsulated=0
  loadct,4
  !p.multi=[0,1,2]
  !P.FONT = 0
  device,/helv,/isolatin1      
  device,font_size=12,/inches,xsize=7.5,ysize=9,xoffset=.5,yoffset=1
  device,filename=iqsweepfit

  ; *** start IQ-loop fitting for both channel 1 (first half of loopdata) and channel 2 (second half of loopdata)
  ; *** this calls a second IDL program ('resfit.pro') to fit the IQ-loop, 'resfit.pro' has to be located in the same subdirectory as this program to be found by IDL
  ; *** resfit.pro saves its results in 2 file, one of it ('series-ps.dat') will be read in later to get the center of the I/Q-loop in order to
  ; *** be able to calculate phase and amplitude from I and Q, ResFit doesn't return any parameters.
  r1 = ResFit(loopdata[*,0:m], Iz1, Qz1,(ts+te)/2.0, StartRes*2-1, StartAtten, outpath, data1a[*,0:m], Iz1a, Qz1a, bg)
  r2 = ResFit(loopdata[*,m+1:lines-1], Iz2, Qz2, (ts+te)/2.0, StartRes*2, StartAtten, outpath, data1a[*,m+1:lines-1], Iz2a, Qz2a, bg)

  ; noise analysis, not used at the moment
  ; noisefile = strcompress(datapath+string(fix(StartT)) +'-'+ string(fix(StartRes)) +'a-'+ string(fix(StartAtten)) +'.ns',/remove_all) 
  ; s1 = NoiseAna(0, noisefile, 0, r1, loopdata[*,0:m], loss)
  ; s2 = NoiseAna(1, noisefile, 1, r2, loopdata[*,m+1:lines-1], loss)     

  device,/close  ; closes PS-output file
endif else print,'IQ-loop fit already present'

; start the fitting procedure 'FitPulses' to fit a generic double exp. decay to all pulses in order to determine the pulse rise time
print,'fitting pulses'
FitPulses,datapath,outpath,pulsename,StartT,StartAtten,StartRes,Npts,sample
print,'done'

end

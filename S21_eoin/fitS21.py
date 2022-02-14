#note, works for IQ files saved using fsweep macro saved on the VNA
#import usual packages 
import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import argrelmin
from scipy import optimize
from scipy.optimize import curve_fit
import pylab
import matplotlib.pyplot as plt
import array as arr
import csv 
#from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages
from datetime import date 
import os


#gets todays date. This is used in the name of the results files that are saved. 
today = date.today() 
todaystr = today.strftime("%Y_%m_%d") 

#E.17 from w thesis - will be used later to  amplitude data fit curve
def amplitudeequation(f, A1, A2, A3, A4, Qr, fr):
    return A1 + A2*(f - fr) + ((A3 + A4*(f - fr)) / (1 + 4 * Qr**2  * ((f - fr)/fr)**2 ))  

#E.11 from Gao thesis - can be used to fit phase data
def phaseequation(f, theta0, Qr, fr):
    return -theta0 + 2*np.arctan(2*Qr*(1 - (f/fr))) 


#ask user to enter filename. Uncomment these too lines if you want to prompt the user
#print('Please enter name of .csv file: ')
#filename = input()

#file name. Hard codes the filename. Comment this out if prompting user for name
foldername = 'Sweep_Files/2022_01_31/'    
filename = 'sweep_2021_11_10_4572530000.csv'
fullname = foldername + filename

#start reader and calculate length of file
with open(fullname) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    row_count = sum(1 for row in csv_reader)  # fileObject is your csv.reader
    
number_of_values = row_count - 8 #accounts for extra lines at start and end of file

#declare 3 arrays, for frequency and I and Q values
freqs = np.zeros(number_of_values)      
Ivalues = np.zeros(number_of_values)
Qvalues = np.zeros(number_of_values)


#open reader again, this time to read in data 
with open(fullname) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    line_count = 0
    for row in csv_reader:
        if line_count <= 6: #opening lines, before actual data 
            #print(f'{", ".join(row)}')
            line_count += 1

        elif line_count < row_count - 1 : #this is the actual IQ data
            #print(f'f = \t{row[0]} Hz, I = {row[1]}, Q = {row[2]}')  #prints data, for debugging 
            
            freqs[line_count-7] = row[0] #saves data to array from csv file
            Ivalues[line_count-7] = row[1] #saves data to array from csv file
            Qvalues[line_count-7] = row[2] #saves data to array from csv file
            
            line_count += 1
    
    #print(f'Processed {line_count} lines.') #prints number of lines read out. Used for debugging. 

freqsMHz = freqs / 10**6
#calculate amplitude of IQ data 
amplitudes = np.sqrt((Ivalues)**2 + (Qvalues)**2)

#find phases of raw data relative to origin
#have to do this for 4 different quadrants    
phases = np.arctan2(Qvalues, Ivalues)


#magnitudes = amplitudes**2 
#normalizedmagnitudes = magnitudes / np.max(magnitudes)

normalizedamplitudes = amplitudes / np.max(amplitudes) #normalizes magnitude data relative to max value
squaremagnitudes = normalizedamplitudes**2 #calculates square magnitude 

p0 = np.array([1, 1, 1, 1, 10000, np.mean(freqsMHz)]) #initial guesses for fit
#p0 = np.array([0.965, 0.22, -0.752, 0.348, 16111, 6120])

#need to define sigma array to weight certain points heavier 
resonant_index = np.argmin(squaremagnitudes) #index of resonant point in array

weight_amplitude = np.ones(number_of_values) #set relative error everywhere to 1
#weight_amplitude[resonant_index - 100 : resonant_index + 100] = 0.01 #set error around resonant to smaller 


#fits square magnitude data to amplitude equation
#popt is the array of fitted parameters
#pcov is the covariance of popt
popt, pcov = curve_fit(amplitudeequation, freqsMHz, squaremagnitudes, p0, sigma = weight_amplitude, absolute_sigma=False, maxfev=10000)   


Q = popt[4] #total Q is taken from the fit
Q = abs(Q)

fr = popt[5]  #resonant frequency is taken from the fit
#fr_MHz = fr / 10**3 #resonant frequency in MHz

Qi = Q / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 
Qi = abs(Qi)

Qc = Q / (1 - np.min(normalizedamplitudes)) #calulculates coupling Q using equation 29 from Zmuidzinas review
Qc = abs(Qc)


#prints results to console
#print(f'Q = {round(Q)}, Qc = {round(Qc)}, Qi = {round(Qi)}, fr = {round(fr_MHz, 4)} MHz')

squaremagnitudesdb = 10*np.log10(squaremagnitudes) #converts square magnitude data to dbs

#plots square magnitudes verus frequency in dbs 
#plt.figure()
#plt.scatter(freqs, squaremagnitudesdb, marker='.', s=5)
#plt.title('$|S_{21}|^2$ (dB) vs Frequency of Raw Data and Fit')
#plt.xlabel('Frequency')
#plt.ylabel('$|S_{21}|^2$ (dB)')
#plt.grid()

#plots square magnitude fit in decibels over data
#plt.plot(freqs, 10*np.log10(amplitudeequation(freqs, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])), 'r')

#calculates midpt of y-axis in plot. Used for saving plot to file
#midpt = (np.min(squaremagnitudesdb) + np.max(squaremagnitudesdb)) / 2


#boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
#boxstr = f'Q = {round(Q)}\nQc = {round(Qc)}\nQi = {round(Qi)}\nfr = {round(fr_MHz, 4)} MHz' #contents of text box
#plt.text(np.min(freqs), np.min(squaremagnitudesdb), boxstr, fontsize=12, bbox=boxprops)#adds textbox to plot


#make name of pdf file to save plots ot
#savefilename = filename.replace('.csv', '') + '_Fit_' + todaystr + '.pdf'

#creates pdf to save plots
#pp = PdfPages(savefilename)       #commentig this out  while fixing to fit phase too
#plt.savefig(pp, format='pdf')



#next need to convert this to IQ data and plot versus raw data 
squaremagnitudefit = amplitudeequation(freqsMHz, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) #saves square magnitude fit data

fitamplitude = np.abs(np.sqrt(squaremagnitudefit)) #takes sqrt to get amplitude

fitamplitudescaled = fitamplitude * np.max(amplitudes) #unnormalizes data relative to max value
 


#get centre of loop
#find centre of loop
raw_data_xc = (np.max(Ivalues) + np.min(Ivalues)) / 2
raw_data_yc = (np.max(Qvalues) + np.min(Qvalues)) / 2


#Next step, remove cable delay term from data. 
#To do this, normalize by cable delay loop
#For now, using circle with radius 4, centered on origin
#Will have to use actual data, after next cooldown
cabledelayphase = np.zeros(number_of_values)

#first find total angle, between first and last point 
#find initial phase

initial_phase = np.arctan2(Qvalues[0], Ivalues[0])
final_phase = np.arctan2(Qvalues[number_of_values-1] , Ivalues[number_of_values-1])

total_phase = abs(final_phase - initial_phase)

phase_per_freq = total_phase / (number_of_values - 1)

i = 0
for x in cabledelayphase:    
    cabledelayphase[i] = initial_phase - phase_per_freq*i        
    i += 1


#need to find amplitude of cable delay loop - take average of start and end  of loop
#amplitudes = np.zeros(number_of_values)
amplitudes = np.sqrt(Ivalues**2 + Qvalues**2)
    
cableamplitude = np.max(amplitudes) #takes max value 
cabledelayI = cableamplitude*np.cos(cabledelayphase) #finds I component of cable
cabledelayQ = cableamplitude*np.sin(cabledelayphase) #finds Q component of cable
    
#plt.plot(cabledelayI, cabledelayQ, 'r')

#substract cable delay from IQ data to normalize
Inormalized = Ivalues - cabledelayI #chceck if this is correct?
Qnormalized = Qvalues - cabledelayQ

#amplitudes_normalized = np.zeros(number_of_values)

amplitudes_normalized = np.sqrt(Inormalized**2 + Qnormalized**2)

#next have to move loop to be centred on (0, 0)

xc = (np.max(Inormalized) + np.min(Inormalized)) / 2 
yc = (np.max(Qnormalized) + np.min(Qnormalized)) / 2

#find alpha - the argument of zc, the centre of the normalized circle
alpha = np.arctan2(yc, xc)

translatedcircleI = Inormalized - xc #finds I component of cable
translatedcircleQ = Qnormalized - yc #finds Q component of cable
    
amplitudes_translated = np.sqrt(translatedcircleI**2 + translatedcircleQ**2)

#plt.scatter(translatedcircleI, translatedcircleQ, s = 5, color = 'm')

theta = np.arctan2(translatedcircleQ, translatedcircleI)
r_translated = np.sqrt((translatedcircleI)**2 + (translatedcircleQ)**2)

#now rotate everything by alpha 

#need to rotate every point in this circle by -alpha 
theta_rotated = theta - alpha #rotates by alpha
I_rotated = r_translated*np.cos(theta_rotated) #finds I component of cable
Q_rotated = r_translated*np.sin(theta_rotated) #finds Q component of cable
    
#amplitudes_rotated = np.zeros(number_of_values)
amplitudes_rotated = np.sqrt(I_rotated**2 + Q_rotated**2)

    
#plt.scatter(I_rotated, Q_rotated, s = 5, color = 'y')    
thetafinal = np.unwrap(theta_rotated)

#now, fit this to phase equation 
p_phase = np.array([1, Q, fr])

#weight_phase = np.ones(number_of_values) #set relative error everywhere to 1
#weight_phase[resonant_index - 250: resonant_index + 250] = 0.05 #set error around resonant to smaller 

#popt_phase, pcov_phase = curve_fit(phaseequation, freqsMHz, thetafinal, p_phase, sigma = weight_phase, absolute_sigma=False)   
popt_phase, pcov_phase = curve_fit(phaseequation, freqsMHz, thetafinal, p_phase)

#print(f'theta0 = {popt[0]}, Qr = {popt[1]}, fr = {popt[2]})')

theta0_phase = popt_phase[0]
Qr_phase = popt_phase[1]
fr_phase = popt_phase[2]

#fr_MHz_phase = fr_phase / 10**3

#plt.plot(freqs, phaseequation(freqs, popt_phase[0], popt_phase[1], popt_phase[2]), 'r')        

Q_phase = abs(Qr_phase)

Qi_phase = Q_phase / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 
Qi_phase = abs(Qi_phase)

Qc_phase = Q_phase / (1 - np.min(normalizedamplitudes)) #calulculates instrinsic Q using equation 29 from Zmuidzinas review
Qc_phase = abs(Qc_phase)


#print(f'Amplitude fit: Q = {Q}, Qc = {Qc}, Qi = {Qi}, fr = {fr*10**9}') 
#print(f'Phase fit: Q = {Q_phase}, Qc = {Qc_phase}, Qi = {Qi_phase}, fr = {fr_phase}') 


#final step, convert phase and amplitude fit back to IQ data and replot fit
fitted_phase = phaseequation(freqsMHz, popt_phase[0], popt_phase[1], popt_phase[2])

#first step, get I and Q data and rotate and translate back to position 


#fitted_phase_rotated = np.zeros(number_of_values)
fitted_phase_rotated = fitted_phase + alpha

fitI = np.mean(amplitudes_rotated)*np.cos(fitted_phase) #finds I component of fit
fitQ = np.mean(amplitudes_rotated)*np.sin(fitted_phase) #finds Q component of fit
        
fitIrotated = np.mean(amplitudes_rotated)*np.cos(fitted_phase_rotated) #finds I component of fit
fitQrotated = np.mean(amplitudes_rotated)*np.sin(fitted_phase_rotated) #finds Q component of fit
    
fitItranslated = fitIrotated + xc
fitQtranslated = fitQrotated + yc
    
fitIunnormalized = fitItranslated + cabledelayI
fitQunnormalized = fitQtranslated + cabledelayQ    

#get phase of unnormalize data
fitphaseunnormalized = np.arctan2(fitQunnormalized, fitIunnormalized)    
    
#last thing to do, combine fitted amplitude and phase data    
Ifitfinal = fitamplitudescaled*np.cos(fitphaseunnormalized) #finds I component of cable
Qfitfinal = fitamplitudescaled*np.sin(fitphaseunnormalized) #finds Q component of cable


#os.chdir('../')

#make name of pdf file to save plots ot
savefilename = foldername + filename.replace('.csv', '') + '_Fit_' + todaystr + '.pdf'



#creates pdf to save plots
pp = PdfPages(savefilename)       #commentig this out  while fixing to fit phase too

#plots square magnitudes verus frequency in dbs 
plt.figure()
plt.scatter(freqsMHz, squaremagnitudesdb, marker='.', s=5)
plt.title('Normalized $|S_{21}|^2$ (dB) vs Frequency of Raw Data and Fit')
plt.xlabel('Frequency (MHz)')
plt.ylabel('$|S_{21}|^2$ (dB)')
plt.grid()

#plots square magnitude fit in decibels over data
plt.plot(freqsMHz, 10*np.log10(amplitudeequation(freqsMHz, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])), 'r')

#calculates midpt of y-axis in plot. Used for saving plot to file
midpt = (np.min(squaremagnitudesdb) + np.max(squaremagnitudesdb)) / 2

boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
boxstr = f'Q = {round(Q_phase)}\nQc = {round(Qc_phase)}\nQi = {round(Qi_phase)}\nfr = {round(fr_phase, 4)} MHz' #contents of text box
plt.text(np.min(freqsMHz), np.min(squaremagnitudesdb), boxstr, fontsize=12, bbox=boxprops)#adds textbox to plot

plt.savefig(pp, format='pdf')
#plt.savefig(pp, format='pdf')


#now plot unwrapped, rotated phase data 
plt.figure()
plt.scatter(freqsMHz, thetafinal, marker='.', s=5)
plt.title('Phase vs Frequency')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Phase (rads)')
plt.plot(freqsMHz, phaseequation(freqsMHz, popt_phase[0], popt_phase[1], popt_phase[2]), 'r')
plt.grid()
boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
boxstr = f'Q = {round(Q_phase)}\nQc = {round(Qc_phase)}\nQi = {round(Qi_phase)}\nfr = {round(fr_phase, 4)} GHz' #contents of text box
plt.text(np.min(freqsMHz), np.min(thetafinal), boxstr, fontsize=12, bbox=boxprops)#adds textbox to plot

plt.savefig(pp, format='pdf')
        
#now plot just origianl IQ data and fitted data 
plt.figure()
plt.scatter(Ivalues, Qvalues, marker='.', s=5)
plt.plot(Ifitfinal, Qfitfinal, 'r')
plt.xlabel('I')
plt.ylabel('Q')
plt.title('Raw and Fitted IQ Data')
plt.grid()
#boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
#boxstr = f'Q = {round(Q_phase)}\nQc = {round(Qc_phase)}\nQi = {round(Qi_phase)}\nfr = {round(fr_MHz_phase, 4)} GHz' #contents of text box
#plt.text(np.min(Ifitfinal), np.min(Qfitfinal), boxstr, fontsize=12, bbox=boxprops)#adds textbox to plot

#creates pdf to save plots
#pp = PdfPages(savefilename)       #commentig this out  while fixing to fit phase too
plt.savefig(pp, format='pdf')
pp.close()






















# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:55:46 2019

@author: Eoin
"""

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
from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages
from datetime import date 

#gets todays date. This is used in the name of the results files that are saved. 
today = date.today() 
todaystr = today.strftime("%Y_%m_%d") 

#E.17 from Gao thesis - will be used later to  amplitude data fit curve
def amplitudeequation(f, A1, A2, A3, A4, Qr, fr):
    return A1 + A2*(f - fr) + ((A3 + A4*(f - fr)) / (1 + 4 * Qr**2  * ((f - fr)/fr)**2 ))  

#E.11 from Gao thesis - can be used to fit phase data (not currently using in code)
def phaseequation(f, theta0, Qr, fr):
    return -theta0 + 2*np.arctan(2*Qr*(1 - (f/fr))) 


#ask user to enter filename. Uncomment these too lines if you want to prompt the user
#print('Please enter name of .csv file: ')
#filename = input()

#file name. Hard codes the filename. Comment this out if prompting user for name
filename = 'fsweepIQ3.csv'


#start reader and calculate length of file
with open(filename) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    row_count = sum(1 for row in csv_reader)  # fileObject is your csv.reader
    #print(row_count)

row_count_nums = row_count - 8 #extra lines at start and end of file

#declare 3 arrays, for frequency and I and Q values
freqs = np.zeros(row_count_nums)       #note -8 due to first view rows of csv file
Ivalues = np.zeros(row_count_nums)
Qvalues = np.zeros(row_count_nums)


#open reader again, this time to read in data 
with open(filename) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    
    line_count = 0
    for row in csv_reader:
        if line_count <= 6: #opening lines, before actual data 
            print(f'{", ".join(row)}')
            line_count += 1

        elif line_count < row_count - 1 : #this is the actual IQ data
            print(f'f = \t{row[0]} Hz, I = {row[1]}, Q = {row[2]}')  #prints data, for debugging 
            
            freqs[line_count-7] = row[0] #saves data to array from csv file
            Ivalues[line_count-7] = row[1] #saves data to array from csv file
            Qvalues[line_count-7] = row[2] #saves data to array from csv file
            
            line_count += 1
    
    print(f'Processed {line_count} lines.') #prints number of lines read out. Used for debugging. 

plt.scatter(Ivalues, Qvalues, color = 'b') #plots raw IQ data.
plt.title('Resonator IQ Plot')
plt.xlabel('I')
plt.ylabel('Q')
plt.grid()
#plt.scatter(Ioffset, Qoffset, color = 'r') 


#calculate magnitude of IQ data 
amplitudes = np.zeros(row_count_nums)
i = 0    
for x in Ivalues:    
    amplitudes[i] = np.sqrt((Ivalues[i])**2 + (Qvalues[i])**2)    
    i += 1


#find phases of raw data relative to origin
#have to do this for 4 different quadrants
phases = np.zeros(row_count_nums)    
i = 0

for x in phases:
    if ((Ivalues[i]) >= 0 and (Qvalues[i]) >= 0):
        phases[i] = np.arctan(np.abs((Qvalues[i])/ (Ivalues[i])))

    elif ((Ivalues[i]) < 0 and (Qvalues[i]) >= 0):
        phases[i] = np.pi - np.arctan(np.abs((Qvalues[i])/ (Ivalues[i])))

    elif ((Ivalues[i]) < 0 and (Qvalues[i]) < 0):
        phases[i] = np.pi + np.arctan(np.abs((Qvalues[i])/ (Ivalues[i])))

    elif ((Ivalues[i]) >= 0 and (Qvalues[i]) < 0):
        phases[i] = 2*np.pi - np.arctan(np.abs((Qvalues[i])/ (Ivalues[i])))

    i += 1


#magnitudes = amplitudes**2 
#normalizedmagnitudes = magnitudes / np.max(magnitudes)

normalizedamplitudes = amplitudes / np.max(amplitudes) #normalizes magnitude data relative to max value
squaremagnitudes = normalizedamplitudes**2 #calculates square magnitude 


#plots square of magnitudes versus frequency
plt.figure()
plt.scatter(freqs, squaremagnitudes, marker='.', s=8)
plt.title('$|S_{21}|^2$ vs Frequency')
plt.xlabel('Frequency')
plt.ylabel('$|S_{21}|^2$')
plt.grid()

p0 = np.array([1, 1, 1, 1, 10000, np.mean(freqs)]) #initial guesses for fit

#need to define sigma array to weight certain points heavier 
resonant_index = np.argmin(squaremagnitudes) #this is resonant point

weight = np.ones(row_count_nums) #set relative error everywhere to 1
weight[resonant_index - 5 : resonant_index + 5] = 0.1 #set error around resonant to smaller 


#fits square magnitude data to amplitude equation
#popt is the array of fitted parameters
#pcov is the covariance of popt
popt, pcov = curve_fit(amplitudeequation, freqs, squaremagnitudes, p0, sigma = weight, absolute_sigma=False)   

#plots fitted curve over square magnitude data
plt.plot(freqs, amplitudeequation(freqs, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]), 'r')

Q = popt[4] #total Q is taken from the fit

fr = popt[5] / 10**9 #resonant frequency is taken from the fit. Converts to GHz

Qi = Q / np.min(normalizedamplitudes) #calulculates instrinsic Q using equation 28 from Zmuidzinas review 

Qc = Q / (1 - np.min(normalizedamplitudes)) #calulculates instrinsic Q using equation 29 from Zmuidzinas review

#prints results to console
print(f'Q = {round(Q)}, Qc = {round(Qc)}, Qi = {round(Qi)}, fr = {round(fr, 4)} GHz')

squaremagnitudesdb = 10*np.log10(squaremagnitudes) #converts square magnitude data to dbs

#plots square magnitudes verus frequency in dbs 
plt.figure()
plt.scatter(freqs, squaremagnitudesdb, marker='.', s=8)
plt.title('$|S_{21}|^2$ (dB) vs Frequency of Raw Data and Fit')
plt.xlabel('Frequency')
plt.ylabel('$|S_{21}|^2$ (dB)')
plt.grid()

#plots square magnitude fit in decibels over data
plt.plot(freqs, 10*np.log10(amplitudeequation(freqs, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])), 'r')

#calculates midpt of y-axis in plot. Used for saving plot to file
midpt = (np.min(squaremagnitudesdb) + np.max(squaremagnitudesdb)) / 2


boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5) #properties of text box on plot
boxstr = f'Q = {round(Q)}\nQc = {round(Qc)}\nQi = {round(Qi)}\nfr = {round(fr, 4)} GHz' #contents of text box
plt.text(np.min(freqs), np.min(squaremagnitudesdb), boxstr, fontsize=12, bbox=boxprops)#adds textbox to plot


#make name of pdf file to save plots ot
savefilename = filename.replace('.csv', '') + '_Fit_' + todaystr + '.pdf'

#creates pdf to save plots
pp = PdfPages(savefilename)
plt.savefig(pp, format='pdf')


#next need to convert this to IQ data and plot versus raw data 
fitI = np.zeros(row_count_nums)
fitQ = np.zeros(row_count_nums)

squaremagnitudefit = amplitudeequation(freqs, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) #saves square magnitude fit data

fitamplitude = np.abs(np.sqrt(squaremagnitudefit)) #takes sqrt to get amplitude

fitamplitudescaled = fitamplitude * np.max(amplitudes) #unnormalizes data relative to max value
 
#fits data to phases equation. Not using now. 
#p0phase = np.array([1, Q, fr])
#poptphase, pcovphase = curve_fit(phaseequation, freqs, phases, p0phase, maxfev = 1000)   
#fitphase = phaseequation(freqs, poptphase[0], poptphase[1], poptphase[2])

#plt.plot(freqs, phaseequation(freqs, poptphase[0], poptphase[1], poptphase[2], 'r'))



i = 0
for x in phases :
    
    fitI[i] = fitamplitudescaled[i]*np.cos(phases[i]) #finds I component of fit
    fitQ[i] = fitamplitudescaled[i]*np.sin(phases[i]) #finds Q component of fit
    i += 1

#plots fitted data on IQ plane
plt.figure()
plt.plot(fitI, fitQ, 'r')
plt.scatter(Ivalues, Qvalues, marker='.', s=8 )
plt.title('IQ Plot of Raw Data and Fit')
plt.xlabel('I')
plt.ylabel('Q')
plt.grid()


#pdf = FPDF()
#pdf.add_page()
#pdf.set_font('Arial', 'B', 16)
#pdf.cell(40, 10, 'Hello World!')
#pdf.output('tuto1.pdf', 'F')


plt.savefig(pp, format='pdf')
pp.close()
#closes pdf
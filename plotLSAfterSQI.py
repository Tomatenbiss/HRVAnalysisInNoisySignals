#http://docs.astropy.org/en/stable/stats/lombscargle.html
import sys
from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle

if len(sys.argv) != 3:
	print("Bitte Datei mit dem zu analysierenden Signal angeben!")
else:
	f = open(sys.argv[1], "r")
	th = float(sys.argv[2])
	

dataEcg = []
for line in f:
	lineValues = line.split()
	dataEcg.append(float(lineValues[2]))

#Perform QRS detection
ecgOut = ecg.ecg(signal=dataEcg, sampling_rate=1000., show=False)

#Calculate RR Tachogram
rPeaks = ecgOut[2]
rrTachogram = []
prevPeak = rPeaks[0]
for peak in rPeaks[1:(len(rPeaks))]:
	rrTachogram.append(peak - prevPeak)
	prevPeak = peak

#Calculate median heartbeat template
templatesForCorrCoef = ecgOut[4]
cleanTemplates = templatesForCorrCoef
medianTemplate = [x / len(cleanTemplates) for x in [sum(x) for x in zip(*cleanTemplates)]]

#Calculate correlation coeffcients
corrCoeffs = []
for template in templatesForCorrCoef:
	corrCoeffs.append(pearsonr(template, medianTemplate)[0])


#create RR Tachogram without bad quality heart beats
rrTachogramAfterSqi = []
tPeaks = []
cnt = 1
for peak in rrTachogram:
	if corrCoeffs[cnt] >= th:
		if corrCoeffs[cnt - 1] >= th:
			rrTachogramAfterSqi.append(peak)
			tPeaks.append(float(rPeaks[cnt] / 1000))
	cnt = cnt + 1;

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma








#freq = np.linspace(0, .4, 100)
#power = LombScargle(tPeaks, rrTachogramAfterSqi).power(freq)
#freq, power = LombScargle(tPeaks, rrTachogramAfterSqi).autopower(minimum_frequency=0, maximum_frequency=.4, samples_per_peak=10)
ls = LombScargle(tPeaks, rrTachogramAfterSqi)
freq, power = ls.autopower(minimum_frequency=0, maximum_frequency=.4, normalization='psd')

freqVlf, powerVlf = ls.autopower(minimum_frequency=0, maximum_frequency=.04, normalization='psd')
freqLf, powerLf = ls.autopower(minimum_frequency=0.04, maximum_frequency=.15, normalization='psd')
freqHf, powerHf = ls.autopower(minimum_frequency=.15, maximum_frequency=.4, normalization='psd')

vlf = 0
lf = 0
hf = 0
for point in powerVlf[1:len(powerVlf)]:
	if not point == float('nan') and not point == float('inf'):
		vlf += point
for point in powerLf:
	if not point == float('nan') and not point == float('inf'):
		lf += point
for point in powerHf:
	if not point == float('nan') and not point == float('inf'):
		hf += point
print("VLF: " + str(vlf) + "LF: " + str(lf) + ", HF: " + str(hf) + ", RATIO LF/HF: " + str(lf/hf))
#print(power.unit)
#freq = np.linspace(0, .4, len(power))

	




#power = power * 2
#power = power * power
#power = power * 1000



#power = movingaverage(power, 125)
#freq = np.linspace(0, .4, len(power))
plt.plot(freq, power)
plt.show()
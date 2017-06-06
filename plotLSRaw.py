import sys
from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle

if len(sys.argv) != 2:
	print("Bitte Datei mit dem zu analysierenden Signal angeben!")
else:
	f = open(sys.argv[1], "r")

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

freq = np.linspace(0, 0.4, 1000)

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma



rPeaks /= 1000
ls = LombScargle(rPeaks[1:len(rPeaks)], rrTachogram)
freq, power = ls.autopower(minimum_frequency=0, maximum_frequency=.4, normalization='psd')

freqVlf, powerVlf = ls.autopower(minimum_frequency=0, maximum_frequency=.04, normalization='psd')
freqLf, powerLf = ls.autopower(minimum_frequency=0.04, maximum_frequency=.15, normalization='psd')
freqHf, powerHf = ls.autopower(minimum_frequency=.15, maximum_frequency=.4, normalization='psd')


print(powerVlf)
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



#freq = np.linspace(0, 0.4, 10000)
#power = LombScargle(rPeaks[1:len(rPeaks)], rrTachogram).power(freq)
#power = movingaverage(power, 1000)
#freq = np.linspace(0, 0.4, len(power))
plt.plot(freq, power)
plt.show()
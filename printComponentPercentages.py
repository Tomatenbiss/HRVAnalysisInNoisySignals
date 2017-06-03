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
			tPeaks.append(rPeaks[cnt])
	cnt = cnt + 1;


def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

freq = np.linspace(0, 0.4, 10000)
power = LombScargle(tPeaks, rrTachogramAfterSqi).power(freq)
power = movingaverage(power, 1000)
freq = np.linspace(0, 0.4, len(power))
plt.plot(freq, power)
plt.show()
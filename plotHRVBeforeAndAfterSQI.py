from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle
from scipy.signal import spectral
import sys

#Getting the data
if len(sys.argv) != 3:
	print("Bitte Datei mit dem zu analysierenden Signal sowie einen Grenzwert angeben!")
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
	else:
		rrTachogramAfterSqi.append(None)
		tPeaks.append(rPeaks[cnt])
	cnt = cnt + 1;




fig = plt.figure()
axRaw = fig.add_subplot(311)
axRaw.plot(rPeaks[1:len(rPeaks)], rrTachogram)

axAfterSQI = fig.add_subplot(312)
axAfterSQI.plot(tPeaks[3:len(tPeaks)], rrTachogramAfterSqi[3:len(rrTachogramAfterSqi)])


axAfterSQI = fig.add_subplot(313)
axAfterSQI.plot(rPeaks, corrCoeffs)

plt.show()


#plt.plot(tPeaks, rrTachogramAfterSqi)
#plt.show()

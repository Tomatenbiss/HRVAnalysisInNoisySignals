from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle

#Getting the data
f = open("sampleDataEcgEda.txt", "r")
dataEcg = []

for line in f:
	lineValues = line.split(",")
	dataEcg.append(float(lineValues[1]))



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
templates = templatesForCorrCoef[-10:-1]
medianTemplate = [x / len(templates) for x in [sum(x) for x in zip(*templates)]]

#Calculate correlation coeffcients
corrCoeffs = []
for template in templatesForCorrCoef:
	corrCoeffs.append(pearsonr(template, medianTemplate)[0])

#create RR Tachogram without bad quality heart beats
rrTachogramAfterSqi = []
t = []
tPeaks = []
cnt = 0
for peak in rrTachogram:
	if corrCoeffs[cnt] >= 0.9:
		rrTachogramAfterSqi.append(peak)
		t.append(cnt)
		tPeaks.append(rPeaks[cnt])
	cnt = cnt + 1;

#rrTachogramAfterSqi = []
#cnt = 1
#for peak in rrTachogram:
#	if corrCoeffs[cnt] >= 0.9:
#		rrTachogramAfterSqi.append(peak)
#	else:
#		rrTachogramAfterSqi.append(None)
#	cnt = cnt + 1



#frequency, power = LombScargle(tPeaks, rrTachogramAfterSqi).autopower()
#power = np.convolve(power, np.ones(20)/20)
#power = power[0:(len(frequency))]
#plt.plot(frequency, power)
#plt.xlim(0, 0.4)
#plt.ylim(0, 0.01)
#plt.show()
plt.plot(tPeaks, rrTachogramAfterSqi)
plt.show()

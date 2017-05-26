from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle

f = open("sampleDataEcgEda.txt", "r")
dataEcg = []

for line in f:
	lineValues = line.split(",")
	dataEcg.append(float(lineValues[1]))

ecgOut = ecg.ecg(signal=dataEcg, sampling_rate=1000., show=False)
rPeaks = ecgOut[2]
rrTachogram = []
prevPeak = rPeaks[0]
for peak in rPeaks[2:-1]:
	rrTachogram.append(peak - prevPeak)
	prevPeak = peak


templatesForCorrCoef = ecgOut[4]
templates = templatesForCorrCoef[-10:-1]
medianTemplate = [x / len(templates) for x in [sum(x) for x in zip(*templates)]]

corrCoeffs = []
for template in templatesForCorrCoef:
	corrCoeffs.append(pearsonr(template, medianTemplate)[0])

rrTachogramAfterSqi = []
t = []
cnt = 1
for peak in rrTachogram:
	if corrCoeffs[cnt] >= 0.9:
		rrTachogramAfterSqi.append(peak)
		t.append(cnt)
	cnt = cnt + 1;

#rrTachogramAfterSqi = []
#cnt = 1
#for peak in rrTachogram:
#	if corrCoeffs[cnt] >= 0.9:
#		rrTachogramAfterSqi.append(peak)
#	else:
#		rrTachogramAfterSqi.append(None)
#	cnt = cnt + 1



frequency, power = LombScargle(t, rrTachogramAfterSqi).autopower()
power = np.convolve(power, np.ones(20)/20)
power = power[0:(len(frequency))]
plt.plot(frequency, power)
plt.xlim(0, 0.4)
plt.ylim(0, 0.01)
plt.show()
#plt.plot(rrTachogramAfterSqi)
#plt.show()

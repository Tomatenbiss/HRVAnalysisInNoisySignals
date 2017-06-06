import sys
from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle
import math

f = open(sys.argv[1], "r")
th = float(sys.argv[2])

def getEcgDataFromFile(file=f, delimiter="", positionInCsvFile=2):
	dataEcg = []
	for line in f:
		lineValues = line.split()
		dataEcg.append(float(lineValues[2]))
	return dataEcg

def getEcgSignal(file=f, delimiter="", positionInCsvFile=2):
	return ecg.ecg(signal=getEcgDataFromFile(), sampling_rate=1000., show=False)

def plotEcgSignal(file=f, delimiter="", positionInCsvFile=2):
	return ecg.ecg(signal=getEcgDataFromFile(), sampling_rate=1000., show=True)

def getRRTachogram(rPeaks):
	rrTachogram = []
	prevPeak = rPeaks[0]
	for peak in rPeaks[1:(len(rPeaks))]:
		rrTachogram.append(peak - prevPeak)
		prevPeak = peak
	return (rPeaks, rrTachogram)

def getMedianHeartbeatTemplate(templatesForCorrCoef, start=0, end=-1):
	cleanTemplates = templatesForCorrCoef[start:end]
	medianTemplate = [x / len(cleanTemplates) for x in [sum(x) for x in zip(*cleanTemplates)]]
	return medianTemplate

def getCorrelationCoefficients(templatesForCorrCoef, medianTemplate):
	corrCoeffs = []
	for template in templatesForCorrCoef:
		corrCoeffs.append(pearsonr(template, medianTemplate)[0])
	return corrCoeffs

def getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks):
	print(rPeaks)
	rrTachogramAfterSqi = []
	tPeaks = []
	cnt = 1
	for peak in rrTachogram:
		if corrCoeffs[cnt] >= th:
			if corrCoeffs[cnt - 1] >= th:
				rrTachogramAfterSqi.append(peak)
				tPeaks.append(float(rPeaks[cnt] / 1000))
		cnt = cnt + 1;
	return (tPeaks, rrTachogramAfterSqi)

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi, minFreq=0, maxFreq=.4, norm='psd'):
	ls = LombScargle(tPeaks, rrTachogramAfterSqi)
	freq, power = ls.autopower(minimum_frequency=minFreq, maximum_frequency=maxFreq, normalization=norm)
	return (freq, power)

def getPower(tPeaks, rrTachogramAfterSqi, minFreq=0, maxFreq=.4, norm='psd'):
	component = 0
	power = getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi, minFreq, maxFreq, norm)[1]
	return sumWithNan(power)

def sumWithNan(lst):
	total = 0
	for item in lst:
		if not math.isnan(item):
			total += item
	return total

def getRatioHFLF(tPeaks, rrTachogramAfterSqi):
	lfPower = getPower(tPeaks, rrTachogramAfterSqi, 0.04, 0.15)
	hfPower = getPower(tPeaks, rrTachogramAfterSqi, 0.15, 0.4)
	return (float(hfPower) / float(lfPower))

def plot(x,y):
	plt.plot(x,y)
	plt.show()

def plotLombScarglePeriodogram():
	ecgSignal = getEcgSignal()
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = getCorrelationCoefficients(ecgSignal[4],medianTemplate)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks)
	(freq, power) = getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi)
	print(getRatioHFLF(tPeaks, rrTachogramAfterSqi))
	#plot(freq, power)



plotLombScarglePeriodogram()
#plotEcgSignal()








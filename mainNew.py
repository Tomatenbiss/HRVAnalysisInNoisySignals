import sys
from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle
import math
from random import randint

f = open(sys.argv[1], "r")

def getEcgDataFromFile(file=f, delimiter="", positionInCsvFile=2):
	dataEcg = []
	for line in f:
		lineValues = line.split()
		dataEcg.append(float(lineValues[2]))
	return dataEcg

def getEcgSignal(file=f, delimiter="", positionInCsvFile=2):
	return ecg.ecg(signal=getEcgDataFromFile(), sampling_rate=1000., show=False)

th = float(sys.argv[2])
percenteageForFakeCoeffs = int(sys.argv[3])
ecgSignal = getEcgSignal()




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


def createFakeCorrelationCoefficients(nBeats, percenteage):
	fakeBeats = nBeats * percenteage / 100
	fakeCorrCoeffs = [1] * nBeats
	while checkIfPercenteageReached(fakeCorrCoeffs, percenteage) == False:
		randPos = randint(0, len(fakeCorrCoeffs) - 1)
		fakeCorrCoeffs[randPos] = 0	
	return fakeCorrCoeffs

def checkIfPercenteageReached(lst, percenteage):
	zeros = 0
	for item in lst:
		if item == 0:
			zeros += 1
	if float(zeros) / len(lst) * 100 >= percenteage:
		return True
	else:
		return False

def plotLombScarglePeriodogramWithFakeCorrelationCoefficients(percenteage):	
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = createFakeCorrelationCoefficients(len(ecgSignal[4]), percenteage)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks)
	(freq, power) = getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi)
	return getRatioHFLF(tPeaks, rrTachogramAfterSqi)
	


def getAverageRatioHFLF(nTimes, percenteage):
	cnt = 0
	avgSum = .0
	minVal = 1000.0
	maxVal = .0
	while cnt < nTimes:
		ratio = float(plotLombScarglePeriodogramWithFakeCorrelationCoefficients(percenteageForFakeCoeffs))
		avgSum += ratio
		if ratio < minVal:
			minVal = ratio
		if ratio > maxVal:
			maxVal = ratio
		cnt += 1
	avgSum /= nTimes
	print("Average: " + str(avgSum) + ", MIN: " + str(minVal) + ", MAX: " + str(maxVal))
		


getAverageRatioHFLF(100, percenteageForFakeCoeffs)
#createFakeCorrelationCoefficients(20, 11)
#plotLombScarglePeriodogram()
#plotEcgSignal()








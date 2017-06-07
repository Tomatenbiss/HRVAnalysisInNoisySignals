import sys
from biosppy.signals import ecg
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from astropy.stats import LombScargle
import math
from random import randint

f = open(sys.argv[1], "r")
dataEcg = []
dataX = []

def getEcgDataFromFile(file=f, delimiter="", positionInCsvFile=2):
	for line in f:
		lineValues = line.split()
		dataEcg.append(float(lineValues[2]))
		dataX.append(float(lineValues[0]))
	return dataEcg


dataEcg = getEcgDataFromFile()

def getEcgSignal(file=f, delimiter="", positionInCsvFile=2):
	return ecg.ecg(signal=dataEcg, sampling_rate=1000., show=False)

th = float(sys.argv[2])
if len(sys.argv) >= 4:
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
				tPeaks.append(float(float(rPeaks[cnt]) / 1000))
		cnt = cnt + 1;
	return (tPeaks, rrTachogramAfterSqi)

def getRRTachogramAfterSQIWithNones(rrTachogram, corrCoeffs, rPeaks):
	rrTachogramAfterSqi = []
	tPeaks = []
	cnt = 1
	for peak in rrTachogram:
		if corrCoeffs[cnt] >= th:
			if corrCoeffs[cnt - 1] >= th:
				rrTachogramAfterSqi.append(peak)
				tPeaks.append(float(float(rPeaks[cnt]) / 1000))
			else:
				rrTachogramAfterSqi.append(None)
				tPeaks.append(float(float(rPeaks[cnt]) / 1000))
		else:
			rrTachogramAfterSqi.append(None)
			tPeaks.append(float(float(rPeaks[cnt]) / 1000))
		cnt = cnt + 1
	return (tPeaks, rrTachogramAfterSqi)

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def getLombScarglePeriodogram(peaks, tachogram, minFreq=0, maxFreq=.4, norm='standard'):
	ls = LombScargle(peaks, tachogram)
	freq, power = ls.autopower(minimum_frequency=minFreq, maximum_frequency=maxFreq, normalization=norm)
	print(sumWithNan(power))
	return (freq, power)

def getPower(tPeaks, rrTachogramAfterSqi, minFreq=0, maxFreq=.4, norm='standard'):
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

def getComponent(tPeaks, rrTachogramAfterSqi, lower, upper):
	vlfPower = getPower(tPeaks, rrTachogramAfterSqi, lower, upper)
	return (vlfPower)

def plot(x,y):
	plt.plot(x,y)
	plt.show()

def plotLombScarglePeriodogram():
	#ecgSignal = getEcgSignal()
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = getCorrelationCoefficients(ecgSignal[4],medianTemplate)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks)
	(freq, power) = getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi)

	#print(getRatioHFLF(tPeaks, rrTachogramAfterSqi))
	#power = movingaverage(power, 100)
	#freq = np.linspace(0, .4, len(power))
	plot(freq, power)

def plotLombScarglePeriodogramRaw():
	#ecgSignal = getEcgSignal()
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	rPeaksCorrected = []
	rPeaks = rPeaks.tolist()
	for peak in rPeaks:
		rPeaksCorrected.append(float(float(peak) / 1000))
	(freq, power) = getLombScarglePeriodogram(rPeaksCorrected[1:len(rPeaks)], rrTachogram)
	#print(getRatioHFLF(tPeaks, rrTachogramAfterSqi))
	plot(freq, power)


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

def plotLombScarglePeriodogramWithFakeCorrelationCoefficientsVLFToHFLF(percenteage):	
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = createFakeCorrelationCoefficients(len(ecgSignal[4]), percenteage)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks)
	(freq, power) = getLombScarglePeriodogram(tPeaks, rrTachogramAfterSqi)
	ratioHFLF = getRatioHFLF(tPeaks, rrTachogramAfterSqi)
	vlf = getComponent(tPeaks, rrTachogramAfterSqi, 0.0, 0.04)
	lf = getComponent(tPeaks, rrTachogramAfterSqi, 0.04, 0.15)
	hf = getComponent(tPeaks, rrTachogramAfterSqi, 0.15, 0.4)
	return (vlf / (lf + hf))
	
def getAverageRatioVLFToHFLF(nTimes, percenteage):
	cnt = 0
	lst = []
	while cnt < nTimes:
		ratio = float(plotLombScarglePeriodogramWithFakeCorrelationCoefficientsVLFToHFLF(percenteageForFakeCoeffs))
		lst.append(ratio)
		cnt += 1
	print(np.mean(lst))
	print(np.std(lst))


def getAverageRatioHFLF(nTimes, percenteage):
	cnt = 0
	lst = []
	while cnt < nTimes:
		ratio = float(plotLombScarglePeriodogramWithFakeCorrelationCoefficients(percenteageForFakeCoeffs))
		lst.append(ratio)
		cnt += 1
	print(np.mean(lst))
	print(np.std(lst))
		
def plotRRTachogramAfterSQI():
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = getCorrelationCoefficients(ecgSignal[4],medianTemplate)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQI(rrTachogram, corrCoeffs, rPeaks)
	plot(tPeaks, rrTachogramAfterSqi)

def plotRRTachogram():
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	plot(rPeaks[1:len(rPeaks)], rrTachogram)

def plotRRTachogramAfterSQIWithNones():
	print("WTF")
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = getCorrelationCoefficients(ecgSignal[4],medianTemplate)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQIWithNones(rrTachogram, corrCoeffs, rPeaks)
	plot(tPeaks, rrTachogramAfterSqi)

def getROCValue():
	rPeaks = ecgSignal[2]
	(rPeaks,rrTachogram) = getRRTachogram(rPeaks)
	medianTemplate = getMedianHeartbeatTemplate(ecgSignal[4])
	corrCoeffs = getCorrelationCoefficients(ecgSignal[4],medianTemplate)
	(tPeaks, rrTachogramAfterSqi) = getRRTachogramAfterSQIWithNones(rrTachogram, corrCoeffs, rPeaks)
	f = open("falsePositives", "r")
	annotations = []
	for line in f:
		items = line.split()
		annotations.append(items[1])
	
	
	TP = 0
	FP = 0
	cnt = 0
	for peak in tPeaks:
		sqiResult = rrTachogramAfterSqi[cnt]
		annotation = annotations[cnt]
		if sqiResult != None:
			if annotation == 'P':
				TP += 1
			else:
				FP += 1
		cnt += 1
	print("TP: " + str(TP) + ", FP: " + str(FP))

		
getROCValue()

def plotPlainEcgSignal():
	x = []
	offset = dataX[0]
	for point in dataX:
		x.append(long((float(point) - offset) * 1000))
	plot(x, dataEcg)

#plotPlainEcgSignal()
#plotRRTachogramAfterSQIWithNones()
#plotRRTachogram()

#getAverageRatioVLFToHFLF(1000, percenteageForFakeCoeffs)
#getAverageRatioHFLF(1000, percenteageForFakeCoeffs)
#createFakeCorrelationCoefficients(20, 11)
#plotLombScarglePeriodogram()
#plotLombScarglePeriodogramRaw()
#plotEcgSignal()


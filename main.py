#!/usr/bin/env python3
import csv
from scipy.interpolate import interp1d
import numpy as np
from math import sqrt

# returns the Radius and the known theta Sin pairs (each as separate np.array)
def readCsvFile(filename):

    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        
        # Get the radius
        line = next(csvReader)
        R = float(line[2])

        # Skip second line
        next(csvReader)
        
        # read the rest
        tmpArray = [(float(row[1]), float(row[2])) for row in csvReader]
        thetas, sinValues = zip(*tmpArray)
       
    return R, thetas, sinValues

# return exmple: 270.5 -> 270 deg & 30 min
def DecimalDegreeToIndividualAngleUnits(decimalDeg):
    degrees = int(decimalDeg)
    minutes = (decimalDeg - degrees) * 60
    seconds = (minutes - int(minutes)) * 60

    return degrees, minutes, seconds

# returns the decimal degrees
def IndividualAngleUnitsToDecimalDegree(degrees, minutes, seconds=0):
    tmpMinutes = minutes + seconds / 60.

    return degrees + tmpMinutes / 360.

def getSizeEpicycle(size_at_0, size_at_90, sinTable, r_for_sin, decimalAngle):
    return size_at_0 + (size_at_90 - size_at_0) * sinTable(decimalAngle) / r_for_sin

def getDecimalAngleFromRotation(revolutionSpeed, elapsedDays, period):
    numRevolutions = (revolutionSpeed * elapsedDays) / (1. * period)
    return (numRevolutions - int(numRevolutions)) * 360

##########

# setup the Sin tables
radiusDeferent, thetas, sinValues = readCsvFile('SinTables/AryabhatiyaSinTable.csv')

InterpolatedSinTable = interp1d(thetas, sinValues)
InterpolatedInverseSinTable = interp1d(sinValues, thetas)

yuga = 4320000
days_in_yuga = 1577917828
days_since_epoch = 1870110 # CHANGE for other dates

# ! ALL OF THE FOLLOWING NUMBERS ARE FOR JUPITER ! #

# revolutions are per yuga
meanPlanet_revolutions = 364220
fast_apogee_revolutions = 4320000

longitude_slow_apogee = IndividualAngleUnitsToDecimalDegree(171, 18)

# size is the circumference, when the deferent has a circumference of 360
# the epicycle has also the sames sizes as the deg + 180
sizeSlow_at_0 = 33
sizeSlow_at_90 = 32
sizeFast_at_0 = 70
sizeFast_at_90 = 72

# 4 step procedure, from surysiddhanta

# 0th step
# calculate the mean planets longitude (lambda_bar)
lambda_bar = getDecimalAngleFromRotation(meanPlanet_revolutions, days_since_epoch, days_in_yuga)
print('lambda_bar: ' + str(lambda_bar))

# 1st step
# apply half the fast equation to the mean planet
lambda_sigma = getDecimalAngleFromRotation(fast_apogee_revolutions, days_since_epoch, days_in_yuga)
print('lambda_sigma: ' + str(lambda_sigma))
kappa_sigma_1 = lambda_bar - lambda_sigma

# constrain kappa, so it's usefull for the sintable
kappaSign = 1.

if(kappa_sigma_1 < 0):
    kappa_sigma_1 += 360

if(kappa_sigma_1 > 180):
    kappaSign = -1
    kappa_sigma_1 -= 180

if(kappa_sigma_1 > 90):
    kappa_sigma_1 -= 90

print('kappa_sigma_1: ' + str(kappa_sigma_1))

sizeFast = getSizeEpicycle(sizeFast_at_0, sizeFast_at_90, InterpolatedSinTable, radiusDeferent, kappa_sigma_1)
print('sizeFast: ' + str(sizeFast))
radiusFast = sizeFast / 360. * radiusDeferent
print('radiusFast: ' + str(radiusFast))

# pretty sure till here it's okay, but check beneath

VMA = kappaSign * InterpolatedSinTable(kappa_sigma_1)
print('VMA: ' + str(VMA))
print('OVM: ' + str(radiusDeferent))
print('VMV: ' + str(radiusFast))

VB = (radiusDeferent * VMA) / radiusFast
print('VB: ' + str(VB))
VB2 = VB * VB
inner = radiusDeferent + sqrt(radiusFast**2 - VB**2)
print('VB2: ' + str(VB2))
print('inner2: ' + str(inner * inner))
radialDistance = sqrt(VB2 + inner * inner)
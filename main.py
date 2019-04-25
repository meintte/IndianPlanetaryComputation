#!/usr/bin/env python3
import csv
from scipy.interpolate import interp1d
import numpy as np
from math import sqrt

class AngleAndSignHandler:
    def __init__(self, thetas, sinValues):
        self.InterpolatedSinTable = interp1d(thetas, sinValues)
        self.InterpolatedInverseSinTable = interp1d(sinValues, thetas)

    # return exmple: 270.5 -> 270 deg & 30 min
    def DecimalDegreeToIndividualAngleUnits(self, decimalDeg):
        degrees = int(decimalDeg)
        minutes = (decimalDeg - degrees) * 60
        seconds = (minutes - int(minutes)) * 60

        return degrees, int(minutes), int(seconds)

    # returns the decimal degrees
    def IndividualAngleUnitsToDecimalDegree(self, degrees, minutes, seconds=0):
        tmpMinutes = minutes + seconds / 60.

        return degrees + tmpMinutes / 360.

    def getPositveAngle(self, decimalAngle):
        while decimalAngle < 0:
            decimalAngle += 360
        return decimalAngle

    # positivity is required
    def _getQuadrantOfAngle(self, decimalAngle):
        # the qudrants are 0, 1, 2, 3
        if (decimalAngle <= 90):
            return 0
        elif (decimalAngle <= 180):
            return 1
        elif (decimalAngle <= 270):
            return 2
        else:
            return 3

    def sinOf(self, decimalAngle):
        angleForSin = self.getPositveAngle(decimalAngle)
        quadrant = self._getQuadrantOfAngle(angleForSin)
        if (quadrant <= 1):
            sign = 1
        else:
            sign = -1

        angleForSin = angleForSin - quadrant*90

        return sign * self.InterpolatedSinTable(angleForSin)

    def arcsinOf(self, sinValue):
        if (sinValue < 0):
            return -1 * self.InterpolatedInverseSinTable(-sinValue)
        else:
            return self.InterpolatedInverseSinTable(sinValue)


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

def getSizeEpicycle(size_at_0, size_at_90, r_for_sin, handlerAngleSin, decimalAngle):
    return size_at_0 + (size_at_90 - size_at_0) * handlerAngleSin.sinOf(decimalAngle) / r_for_sin

def getDecimalAngleFromRotation(revolutionSpeed, elapsedDays, period):
    numRevolutions = (revolutionSpeed * elapsedDays) / (1. * period)
    return (numRevolutions - int(numRevolutions)) * 360

def getFastEquation(radiusFast, radiusDeferent, handlerAngleSin, kappa):
    
    sinKappa = handlerAngleSin.sinOf(kappa)

    print('VMA: ' + str(sinKappa))
    print('OVM: ' + str(radiusDeferent))
    print('VMV: ' + str(radiusFast))

    VB = (radiusFast * sinKappa) / radiusDeferent
    print('VB: ' + str(VB))
    radialDistance = sqrt(VB**2 + (radiusDeferent + sqrt(sinKappa**2 - VB**2))**2)
    print('radialDistance: ' + str(radialDistance))

    sigma = handlerAngleSin.arcsinOf(radiusFast * sinKappa / radialDistance)
   
    return sigma

##########

# setup the Sin tables
# it's assumed that the sin table only gives vales for angles in [0,90 deg]
radiusDeferent, thetas, sinValues = readCsvFile('SinTables/AryabhatiyaSinTable.csv')

handler = AngleAndSignHandler(thetas, sinValues)

yuga = 4320000
days_in_yuga = 1577917828
days_since_epoch = 1870110 # CHANGE for other dates

# ! ALL OF THE FOLLOWING NUMBERS ARE FOR JUPITER ! #

# revolutions are per yuga
meanPlanet_revolutions = 364220
fast_apogee_revolutions = 4320000

longitude_slow_apogee = handler.IndividualAngleUnitsToDecimalDegree(171, 18)

# size is the circumference, when the deferent has a circumference of 360
# the epicycle has also the sames sizes as the deg + 180
sizeSlow_at_0 = 33
sizeSlow_at_90 = 32
sizeFast_at_0 = 70
sizeFast_at_90 = 72

# 4 step procedure, from suryasiddhanta

# 0th step
# calculate the mean planets longitude (lambda_bar)
lambda_bar = getDecimalAngleFromRotation(meanPlanet_revolutions, days_since_epoch, days_in_yuga)
print('lambda_bar: ' + str(lambda_bar))

# 1st step
# apply half the fast equation to the mean planet
lambda_sigma = getDecimalAngleFromRotation(fast_apogee_revolutions, days_since_epoch, days_in_yuga)
print('lambda_sigma: ' + str(lambda_sigma))
kappa_sigma_1 = lambda_bar - lambda_sigma

print('kappa_sigma_1: ' + str(kappa_sigma_1))

# get the current size of the epicycle
sizeFast = getSizeEpicycle(sizeFast_at_0, sizeFast_at_90, radiusDeferent, handler, kappa_sigma_1)
print('sizeFast: ' + str(sizeFast))
# convert the size to the radius
radiusFast = sizeFast / 360. * radiusDeferent
print('radiusFast: ' + str(radiusFast))


sigma_1 = getFastEquation(radiusFast, radiusDeferent,handler, kappa_sigma_1)
print('sigma_1: ' + str(sigma_1))

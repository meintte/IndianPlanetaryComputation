#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from scipy.interpolate import interp1d
import numpy as np
from math import sqrt

class AngleAndSinHandler:
    def __init__(self, thetas, sinValues):
        self.InterpolatedSinTable = interp1d(thetas, sinValues)
        self.InterpolatedInverseSinTable = interp1d(sinValues, thetas)

    # return exmple: 270.5 -> 270 deg & 30 min
    def DecimalDegreeToIndividualAngleUnits(self, decimalDeg):
        degrees = int(decimalDeg)
        minutes = (decimalDeg - degrees) * 60
        seconds = (minutes - int(minutes)) * 60

        return degrees, int(minutes), seconds

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

    def printAngle(self, name, angle, inDecimal=True):
        if inDecimal:
            print('{:20}: {}°'.format(name, angle))
        else:
            _deg, _min, _sec = self.DecimalDegreeToIndividualAngleUnits(angle)

            print('{:20}: {}° {}\' {}\'\''.format(name, _deg, _min, _sec))


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
    return size_at_0 + (size_at_90 - size_at_0) * abs(handlerAngleSin.sinOf(decimalAngle)) / r_for_sin

def getDecimalAngleFromRotation(revolutionSpeed, elapsedDays, period):
    numRevolutions = (revolutionSpeed * elapsedDays) / (1. * period)
    return (numRevolutions - int(numRevolutions)) * 360

def getFastEquation(radiusFast, radiusDeferent, handlerAngleSin, kappa):
    
    sinKappa = handlerAngleSin.sinOf(kappa)

    #print('VMA: \t' + str(sinKappa))
    #print('OVM: \t' + str(radiusDeferent))
    #print('VMV: \t' + str(radiusFast))

    VB = (radiusFast * sinKappa) / radiusDeferent
    #print('VB: \t' + str(VB))
    radialDistance = sqrt(VB**2 + (radiusDeferent + sqrt(sinKappa**2 - VB**2))**2)
    #print('radialDistance: \t' + str(radialDistance))

    sigma = handlerAngleSin.arcsinOf(radiusFast * sinKappa / radialDistance)
   
    return sigma

def getSlowEquation(radiusSlow, radiusDeferent, handlerAngleSin, kappa):
    
    sinKappa = handlerAngleSin.sinOf(kappa)

    mu = handlerAngleSin.arcsinOf(radiusSlow * sinKappa / radiusDeferent)

    return mu

##########

# setup the Sin tables
# it's assumed that the sin table only gives vales for angles in [0,90 deg]
radiusDeferent, thetas, sinValues = readCsvFile('SinTables/AryabhatiyaSinTable.csv')

handler = AngleAndSinHandler(thetas, sinValues)

# evidence suggest that values are rounded to the nearest minute
doRounding = False
# print all steps
printAll = True
# print angles in decimalDegree
printDecimalDegree = False

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

### END OF PARAMETER DECLARATION

# 4 step procedure, from suryasiddhanta

# 0th step
# calculate the mean planets longitude (lambda_bar)
lambda_bar = getDecimalAngleFromRotation(meanPlanet_revolutions, days_since_epoch, days_in_yuga)
handler.printAngle('lambda_bar', lambda_bar, printDecimalDegree)

################# START 1st step #################
# apply half the fast equation to the mean planet
lambda_sigma = getDecimalAngleFromRotation(fast_apogee_revolutions, days_since_epoch, days_in_yuga)
handler.printAngle('lambda_sigma', lambda_sigma, printDecimalDegree)

kappa_sigma_1 = lambda_bar - lambda_sigma
handler.printAngle('kappa_sigma_1', kappa_sigma_1, printDecimalDegree)
kappa_sigma_1 = handler.getPositveAngle(kappa_sigma_1)
handler.printAngle('kappa_sigma_1', kappa_sigma_1, printDecimalDegree)

# get the current size of the epicycle
sizeFast = getSizeEpicycle(sizeFast_at_0, sizeFast_at_90, radiusDeferent, handler, kappa_sigma_1)
#print('sizeFast: \t' + str(sizeFast))
# convert the size to the radius
radiusFast = sizeFast / 360. * radiusDeferent
#print('radiusFast: \t' + str(radiusFast))

sigma_1 = getFastEquation(radiusFast, radiusDeferent, handler, kappa_sigma_1)
handler.printAngle('sigma_1', sigma_1, printDecimalDegree)
sigma_1 = handler.getPositveAngle(sigma_1)
handler.printAngle('sigma_1', sigma_1, printDecimalDegree)

# plus or minus ?? took minus for now
lambda_1 = lambda_bar - 0.5 * sigma_1
handler.printAngle('lambda_1', lambda_1, printDecimalDegree)

################# END 1st step #################

################# START 2nd step #################
# apply half the slow equation to the computed result
lambda_mu = longitude_slow_apogee
handler.printAngle('lambda_mu', lambda_mu, printDecimalDegree)

kappa_mu_1 = lambda_1 - lambda_mu
handler.printAngle('kappa_mu_1', kappa_mu_1, printDecimalDegree)
kappa_mu_1 = handler.getPositveAngle(kappa_mu_1)
handler.printAngle('kappa_mu_1', kappa_mu_1, printDecimalDegree)

# get the current size of the epicycle
sizeSlow = getSizeEpicycle(sizeSlow_at_0, sizeSlow_at_90, radiusDeferent, handler, kappa_mu_1)
#print('sizeSlow: \t' + str(sizeSlow))
# convert the size to the radius
radiusSlow = sizeSlow / 360. * radiusDeferent
#print('radiusSlow: \t' + str(radiusSlow))

mu_1 = getSlowEquation(radiusSlow, radiusDeferent, handler, kappa_mu_1)
handler.printAngle('mu_1', mu_1, printDecimalDegree)
mu_1 = handler.getPositveAngle(mu_1)
handler.printAngle('mu_1', mu_1, printDecimalDegree)

# plus or minus ?? took plus for now
lambda_2 = lambda_1 + 0.5 * mu_1
handler.printAngle('lambda_2', lambda_2, printDecimalDegree)

################# END 2nd step #################

################# START 3rd step #################
# start form the computed result, compute the slow equation,
# apply it whole to the mean planet
kappa_mu_2 = lambda_2 - lambda_mu
handler.printAngle('kappa_mu_2', kappa_mu_2, printDecimalDegree)
kappa_mu_2 = handler.getPositveAngle(kappa_mu_2)
handler.printAngle('kappa_mu_2', kappa_mu_2, printDecimalDegree)

# get the current size of the epicycle
sizeSlow = getSizeEpicycle(sizeSlow_at_0, sizeSlow_at_90, radiusDeferent, handler, kappa_mu_2)
#print('sizeSlow: \t' + str(sizeSlow))
# convert the size to the radius
radiusSlow = sizeSlow / 360. * radiusDeferent
#print('radiusSlow: \t' + str(radiusSlow))

mu_2 = getSlowEquation(radiusSlow, radiusDeferent, handler, kappa_mu_2)
handler.printAngle('mu_2', mu_2, printDecimalDegree)
mu_2 = handler.getPositveAngle(mu_2)
handler.printAngle('mu_2', mu_2, printDecimalDegree)

# plus or minus ?? took plus for now
lambda_3 = lambda_bar + mu_2
handler.printAngle('lambda_3', lambda_3, printDecimalDegree)

################# END 3rd step #################

################# START 4th step #################
# apply the whole fast equation to the computed result
kappa_sigma_2 = lambda_3 - lambda_sigma
handler.printAngle('kappa_sigma_2', kappa_sigma_2, printDecimalDegree)
kappa_sigma_2 = handler.getPositveAngle(kappa_sigma_2)
handler.printAngle('kappa_sigma_2', kappa_sigma_2, printDecimalDegree)

# get the current size of the epicycle
sizeFast = getSizeEpicycle(sizeFast_at_0, sizeFast_at_90, radiusDeferent, handler, kappa_sigma_2)
#print('sizeFast: \t' + str(sizeFast))
# convert the size to the radius
radiusFast = sizeFast / 360. * radiusDeferent
#print('radiusFast: \t' + str(radiusFast))

sigma_2 = getFastEquation(radiusFast, radiusDeferent, handler, kappa_sigma_2)
handler.printAngle('sigma_2', sigma_2, printDecimalDegree)
sigma_2 = handler.getPositveAngle(sigma_2)
handler.printAngle('sigma_2', sigma_2, printDecimalDegree)

# plus or minus ?? took minus for now
lambda_true = lambda_3 - sigma_2
handler.printAngle('lambda_true', lambda_true, printDecimalDegree)
lambda_true = handler.getPositveAngle(lambda_true)
handler.printAngle('sigma_2', lambda_true, printDecimalDegree)

################# END 4th step #################
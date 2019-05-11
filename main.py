#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from math import sqrt

import numpy as np
import yaml
from scipy.interpolate import interp1d


class PlanetVariables:
    def __init__(self, name):
        self.name = name
        self._getDataFromFile()

    def _getDataFromFile(self):
        with open("planetVariables/" + self.name + ".yml", 'r') as ymlfile:
            yamlFile = yaml.load(ymlfile)

        self.meanPlanet_revolutions = yamlFile['mean planet revolutions per yuga']

        self.longitude_slow_apogee = yamlFile['longitude slow apogee']
        self.sizeSlow_at_0 = yamlFile['size slow epicycle at 0 & 180']
        self.sizeSlow_at_90 = yamlFile['size slow epicycle at 90 & 270']

        if self.name != 'Sun':
            self.fast_apogee_revolutions = yamlFile['fast apogee revolutions per yuga']
            self.sizeFast_at_0 = yamlFile['size fast epicycle at 0 & 180']
            self.sizeFast_at_90 = yamlFile['size fast epicycle at 90 & 270']


class AngleAndSinHandler:
    def __init__(self, thetas, sinValues):
        self.InterpolatedSinTable = interp1d(thetas, sinValues)
        self.InterpolatedInverseSinTable = interp1d(sinValues, thetas)

    # return exmple: 270.5 -> 270 deg & 30 min
    # bug with input 30.2666667, the int in min makes from 16.0 15...
    def DecimalDegreeToIndividualAngleUnits(self, decimalDeg):
        degrees = int(decimalDeg)
        minutes = int((decimalDeg - degrees) * 60)
        seconds = (decimalDeg - degrees - minutes/60.)*3600

        return degrees, minutes, seconds

    # returns the decimal degrees
    def IndividualAngleUnitsToDecimalDegree(self, degrees, minutes, seconds=0):
        tmpMinutes = minutes + seconds / 60.

        return degrees + tmpMinutes / 60.

    def getPositveAngle(self, decimalAngle):
        while decimalAngle < 0:
            decimalAngle += 360
        while decimalAngle > 360:
            decimalAngle -= 360
        return decimalAngle

    def roundToMinutes(self, decimalAngle):
        _deg, _min, _sec = self.DecimalDegreeToIndividualAngleUnits(
            decimalAngle)
        if (_sec >= 30.):
            return self.IndividualAngleUnitsToDecimalDegree(_deg, _min+1, 0)
        else:
            return self.IndividualAngleUnitsToDecimalDegree(_deg, _min, 0)

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
        # the quadrant numberation goes from 0 to 3
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

    def printAngle(self, name, decimalAngle, inDecimal=True):
        if inDecimal:
            print('{:20}: {}°'.format(name, decimalAngle))
        else:
            _deg, _min, _sec = self.DecimalDegreeToIndividualAngleUnits(
                decimalAngle)

            print('{:20}: {}° {}\' {}\'\''.format(name, _deg, _min, _sec))

    def makePositiveRoundAndPrint(self, name, angle, inDecimal=True, doRound=False):
        # make positive
        angle = self.getPositveAngle(angle)
        # do the rounding
        if doRound:
            angle = self.roundToMinutes(angle)
        # print the angle
        self.printAngle(name, angle, inDecimal)

        return angle


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
    return size_at_0 + (size_at_90 - size_at_0) * abs(handlerAngleSin.sinOf(decimalAngle)) / (1. * r_for_sin)


def getRadiusEpicycle(size_at_0, size_at_90, radiusDeferent, handlerAngleSin, decimalAngle, printAll):
    sizeEpicycle = getSizeEpicycle(
        size_at_0, size_at_90, radiusDeferent, handlerAngleSin, decimalAngle)

    radiusEpicycle = sizeEpicycle / 360. * radiusDeferent

    if printAll:
        print('{:20}: {}'.format('sizeEpicycle', sizeEpicycle))
        print('{:20}: {}'.format('radiusEpicycle', radiusEpicycle))

    return radiusEpicycle


def getDecimalAngleFromRotation(revolutionSpeed, elapsedDays, period):
    numRevolutions = (revolutionSpeed * elapsedDays) / (1. * period)
    return (numRevolutions - int(numRevolutions)) * 360


def getFastEquation(radiusFast, radiusDeferent, handlerAngleSin, kappa, printAll):

    sinKappa = handlerAngleSin.sinOf(kappa)

    VB = (radiusFast * sinKappa) / radiusDeferent

    radialDistance = sqrt(
        VB**2 + (radiusDeferent + sqrt(sinKappa**2 - VB**2))**2)

    sigma = handlerAngleSin.arcsinOf(radiusFast * sinKappa / radialDistance)

    if printAll:
        print('{:20}: {}'.format('sinKappa', sinKappa))
        print('{:20}: {}'.format('radiusDeferent', radiusDeferent))
        print('{:20}: {}'.format('radiusFast', radiusFast))
        print('{:20}: {}'.format('radialDistance', radialDistance))

    return sigma


def getSlowEquation(radiusSlow, radiusDeferent, handlerAngleSin, kappa, printAll):

    sinKappa = handlerAngleSin.sinOf(kappa)

    if printAll:
        print('{:20}: {}'.format('sinKappa', sinKappa))

    mu = handlerAngleSin.arcsinOf(radiusSlow * sinKappa / radiusDeferent)

    return mu

#############################################################


def doSunProcedure(
        _yuga, _days_in_yuga, _days_since_epoch,
        _radiusDeferent, _handler,
        _meanPlanet_revolutions, _longitude_slow_apogee,
        _sizeSlow_at_0, _sizeSlow_at_90,
        _doRounding, _printDecimalDegree, _printAll):

    # mean planet calculation
    lambda_bar = getDecimalAngleFromRotation(
        _meanPlanet_revolutions, _days_since_epoch, _days_in_yuga)
    lambda_bar = _handler.makePositiveRoundAndPrint(
        'lambda_bar', lambda_bar, _printDecimalDegree, _doRounding)

    # apply half the slow equation to the computed result
    lambda_mu = _longitude_slow_apogee
    lambda_mu = _handler.makePositiveRoundAndPrint(
        'lambda_mu', lambda_mu, _printDecimalDegree, _doRounding)

    kappa_mu = lambda_bar - lambda_mu
    kappa_mu = _handler.makePositiveRoundAndPrint(
        'kappa_mu', kappa_mu, _printDecimalDegree, _doRounding)

    # get the current radius of the epicycle
    radiusSlow = getRadiusEpicycle(
        _sizeSlow_at_0, _sizeSlow_at_90, _radiusDeferent, _handler, kappa_mu, _printAll)

    mu = getSlowEquation(radiusSlow, _radiusDeferent,
                         _handler, kappa_mu, _printAll)
    mu = _handler.makePositiveRoundAndPrint(
        'mu', mu, _printDecimalDegree, _doRounding)

    # plus or minus? use the secondSign...
    lambda_true = lambda_bar + mu
    lambda_true = _handler.makePositiveRoundAndPrint(
        'lambda_true', lambda_true, _printDecimalDegree, _doRounding)

#############################################################


def do4stepProcedure(
        _yuga, _days_in_yuga, _days_since_epoch,
        _radiusDeferent, _handler,
        _meanPlanet_revolutions, _fast_apogee_revolutions, _longitude_slow_apogee,
        _sizeSlow_at_0, _sizeSlow_at_90, _sizeFast_at_0, _sizeFast_at_90,
        _doRounding, _printDecimalDegree, _printAll,
        _firstSign, _secondSign, _thirdSign, _fourthSign):

    # 4 step procedure, from suryasiddhanta

    # 0th step
    # calculate the mean planets longitude (lambda_bar)
    lambda_bar = getDecimalAngleFromRotation(
        _meanPlanet_revolutions, _days_since_epoch, _days_in_yuga)
    lambda_bar = _handler.makePositiveRoundAndPrint(
        'lambda_bar', lambda_bar, _printDecimalDegree, _doRounding)

    ################# START 1st step #################
    # apply half the fast equation to the mean planet
    lambda_sigma = getDecimalAngleFromRotation(
        _fast_apogee_revolutions, _days_since_epoch, _days_in_yuga)
    lambda_sigma = _handler.makePositiveRoundAndPrint(
        'lambda_sigma', lambda_sigma, _printDecimalDegree, _doRounding)

    kappa_sigma_1 = lambda_bar - lambda_sigma
    kappa_sigma_1 = _handler.makePositiveRoundAndPrint(
        'kappa_sigma_1', kappa_sigma_1, _printDecimalDegree, _doRounding)

    # get the current radius of the epicycle
    radiusFast = getRadiusEpicycle(
        _sizeFast_at_0, _sizeFast_at_90, _radiusDeferent, _handler, kappa_sigma_1, _printAll)

    sigma_1 = getFastEquation(
        radiusFast, _radiusDeferent, _handler, kappa_sigma_1, _printAll)
    sigma_1 = _handler.makePositiveRoundAndPrint(
        'sigma_1', sigma_1, _printDecimalDegree, _doRounding)

    # plus or minus? use the firstSign...
    lambda_1 = lambda_bar + _firstSign * 0.5 * sigma_1
    lambda_1 = _handler.makePositiveRoundAndPrint(
        'lambda_1', lambda_1, _printDecimalDegree, _doRounding)

    ################# END 1st step #################

    ################# START 2nd step #################
    # apply half the slow equation to the computed result
    lambda_mu = _longitude_slow_apogee
    lambda_mu = _handler.makePositiveRoundAndPrint(
        'lambda_mu', lambda_mu, _printDecimalDegree, _doRounding)

    kappa_mu_1 = lambda_1 - lambda_mu
    kappa_mu_1 = _handler.makePositiveRoundAndPrint(
        'kappa_mu_1', kappa_mu_1, _printDecimalDegree, _doRounding)

    # get the current radius of the epicycle
    radiusSlow = getRadiusEpicycle(
        _sizeSlow_at_0, _sizeSlow_at_90, _radiusDeferent, _handler, kappa_mu_1, _printAll)

    mu_1 = getSlowEquation(radiusSlow, _radiusDeferent,
                           _handler, kappa_mu_1, _printAll)
    mu_1 = _handler.makePositiveRoundAndPrint(
        'mu_1', mu_1, _printDecimalDegree, _doRounding)

    # plus or minus? use the secondSign...
    lambda_2 = lambda_1 + _secondSign * 0.5 * mu_1
    lambda_2 = _handler.makePositiveRoundAndPrint(
        'lambda_2', lambda_2, _printDecimalDegree, _doRounding)

    ################# END 2nd step #################

    ################# START 3rd step #################
    # start form the computed result, compute the slow equation,
    # apply it whole to the mean planet
    kappa_mu_2 = lambda_2 - lambda_mu
    kappa_mu_2 = _handler.makePositiveRoundAndPrint(
        'kappa_mu_2', kappa_mu_2, _printDecimalDegree, _doRounding)

    # get the current radius of the epicycle
    radiusSlow = getRadiusEpicycle(
        _sizeSlow_at_0, _sizeSlow_at_90, _radiusDeferent, _handler, kappa_mu_2, _printAll)

    mu_2 = getSlowEquation(radiusSlow, _radiusDeferent,
                           _handler, kappa_mu_2, _printAll)
    mu_2 = _handler.makePositiveRoundAndPrint(
        'mu_2', mu_2, _printDecimalDegree, _doRounding)

    # plus or minus? use the thridSign...
    lambda_3 = lambda_bar + _thirdSign * mu_2
    lambda_3 = _handler.makePositiveRoundAndPrint(
        'lambda_3', lambda_3, _printDecimalDegree, _doRounding)

    ################# END 3rd step #################

    ################# START 4th step #################
    # apply the whole fast equation to the computed result
    kappa_sigma_2 = lambda_3 - lambda_sigma
    kappa_sigma_2 = _handler.makePositiveRoundAndPrint(
        'kappa_sigma_2', kappa_sigma_2, _printDecimalDegree, _doRounding)

    # get the current size of the epicycle
    radiusFast = getRadiusEpicycle(
        _sizeFast_at_0, _sizeFast_at_90, _radiusDeferent, _handler, kappa_sigma_2, _printAll)

    sigma_2 = getFastEquation(
        radiusFast, _radiusDeferent, _handler, kappa_sigma_2, _printAll)
    sigma_2 = _handler.makePositiveRoundAndPrint(
        'sigma_2', sigma_2, _printDecimalDegree, _doRounding)

    # plus or minus? use the fourthSign...
    lambda_true = lambda_3 + _fourthSign * sigma_2
    lambda_true = _handler.makePositiveRoundAndPrint(
        'lambda_true', lambda_true, _printDecimalDegree, _doRounding)

    ################# END 4th step #################

#############################################################


def allPosibilityWay(yuga, days_in_yuga, days_since_epoch, radiusDeferent, handler, planet, doRounding, printDecimalDegree, printAll):
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                for l in [-1, 1]:
                    print(
                        "####################################### " '{},{},{},{}'.format(i, j, k, l))
                    do4stepProcedure(
                        yuga, days_in_yuga, days_since_epoch,
                        radiusDeferent, handler,
                        planet.meanPlanet_revolutions, planet.fast_apogee_revolutions, planet.longitude_slow_apogee,
                        planet.sizeSlow_at_0, planet.sizeSlow_at_90, planet.sizeFast_at_0, planet.sizeFast_at_90,
                        doRounding, printDecimalDegree, printAll,
                        i, j, k, l)

#############################################################


if __name__ == "__main__":

    # get the global variables from the yaml file
    with open("globalVariables.yml", 'r') as ymlfile:
        globalVars = yaml.load(ymlfile)

    # setup the Sin tables
    # it's assumed that the sin table only gives vales for angles in [0,90 deg]
    radiusDeferent, thetas, sinValues = readCsvFile(
        globalVars['sin table'])

    handler = AngleAndSinHandler(thetas, sinValues)

    # evidence suggest that angle values are rounded to the nearest minute
    doRounding = globalVars['round to minutes']

    # print angles in decimalDegree
    printDecimalDegree = globalVars['print in decimal degrees']

    # print all steps
    printAll = globalVars['print all steps']

    yuga = globalVars['yuga']
    days_in_yuga = globalVars['days in a yuga']
    days_since_epoch = globalVars['days since the epoch']

    planets = []

    for planetToCalculate in globalVars['do calculations for']:
        if planetToCalculate == 'Sun':
            planets.append(PlanetVariables(planetToCalculate))
        elif planetToCalculate == 'Moon':
            print(planetToCalculate + "is not yet implemented...")
        elif planetToCalculate == 'Mars' or planetToCalculate == 'Mercury' or planetToCalculate == 'Jupiter' or planetToCalculate == 'Venus' or planetToCalculate == 'Saturn':
            planets.append(PlanetVariables(planetToCalculate))
        else:
            print("Unknown planet! Please check for typos")

    for p in planets:

        print(p.name)
        print("")

        if(p.name == 'Sun'):
            doSunProcedure(
                yuga, days_in_yuga, days_since_epoch,
                radiusDeferent, handler,
                p.meanPlanet_revolutions, p.longitude_slow_apogee,
                p.sizeSlow_at_0, p.sizeSlow_at_90,
                doRounding, printDecimalDegree, printAll)
        else:
            do4stepProcedure(
                yuga, days_in_yuga, days_since_epoch,
                radiusDeferent, handler,
                p.meanPlanet_revolutions, p.fast_apogee_revolutions, p.longitude_slow_apogee,
                p.sizeSlow_at_0, p.sizeSlow_at_90, p.sizeFast_at_0, p.sizeFast_at_90,
                doRounding, printDecimalDegree, printAll,
                -1, 1, 1, -1)

        print("")

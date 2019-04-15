#!/usr/bin/env python3
import csv

# returns the Radius, SinTable and dSin from a File
def readCsvFile(filename):

    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        
        # Get the radius:
        line = next(csvReader)
        R = line[2]

        # Skip second line:
        next(csvReader)
        
        # read the rest
        sinTable = {rows[1]:rows[2] for rows in csvReader}
        dsin = {rows[1]:rows[3] for rows in csvReader}
 
    return R, sinTable, dsin
 
 
Radius, SinTable, DeltaSin = readCsvFile('SinTables/AryabhatiyaSinTable.csv')

print("SinTable for R = " + Radius)

for x, sin in SinTable.items():
    print('{}:\t{}'.format(x, sin))
import csv
import numpy as np

def saveDarkTimesCSV(objectIDs,dts,filename):
    with open( filename, 'w') as outfile:
        writer = csv.writer( outfile )
        for id, darktimes in zip(objectIDs,dts):
            writer.writerow([id]+darktimes.tolist())

def readDarkTimesCSV(filename):
    with open( filename, 'r') as infile:
        reader = csv.reader( infile )
        objectIDs = []
        dts = []
        for row in reader:
            objectIDs.append(int(row[0]))
            dt = np.array(row[1:],dtype='f')
            dts.append(dt)
    return (objectIDs,dts)

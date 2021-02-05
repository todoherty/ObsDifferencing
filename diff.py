#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import sys
import os
from datetime import datetime
from casacore.tables import table, taql, maketabdesc, makescacoldesc
from argparse import ArgumentParser
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.utils import data as astropydata 
import astropy.units as u

#astropydata.import_download_cache('astropy_data_2020-08-25.cache')

#Root directory
print(sys.argv)
#root = "/astro/mwasci/todoherty/"
root = sys.argv[3]

#print(EarthLocation.get_site_names())
#import astropy
#astropy.coordinates.earth.EarthLocation._get_site_registry(force_download=True)
#print(EarthLocation.get_site_names())

#Opening ms files
print(sys.argv)
obsid1 = sys.argv[1]
obsid2 = sys.argv[2]

obs1 = root + sys.argv[1]
obs2 = root + sys.argv[2]

mset1 = table("{0}/{1}.ms".format(obs1,obsid1), readonly=True)
mset2 = table("{0}/{1}.ms".format(obs2,obsid2), readonly=True)



#ants to flag
ants1 = set(mset1.getcol("ANTENNA1"))
ants2 = set(mset2.getcol("ANTENNA1"))

ants_to_flag = np.array(list(ants1.symmetric_difference(ants2)))

if(ants_to_flag.size != 0):
	print("Antennas to flag: " + str(ants_to_flag))
else:
	print("No antennas to flag")


#EarthLocation._get_site_registry(force_download=True)

#time to flag
time1 = list(set(mset1.getcol("TIME")))
time2 = list(set(mset2.getcol("TIME")))

print(time1)
t1_mjd = [ x/(24*3600) for x in time1]
t2_mjd = [ x/(24*3600) for x in time2]
print(t1_mjd)


#mwa = EarthLocation.of_site("Murchison Widefield Array")
mwa = EarthLocation.from_geodetic(lon=116.67081524, lat=-26.7033194, height=377.83, ellipsoid='WGS84')
print(mwa.geocentric)
print(mwa.geodetic)
time1_mjd = Time(t1_mjd, format='mjd', scale='utc', location=mwa)
time2_mjd = Time(t2_mjd, format='mjd', scale='utc', location=mwa)

print(time1_mjd)

#sidereal_time() returns hrs, *3600 for seconds
lst1 = time1_mjd.sidereal_time('apparent', model = 'IAU2006A') * 3600
lst2 = time2_mjd.sidereal_time('apparent', model = 'IAU2006A') * 3600

#finding any times that need to be flagged from both ms
time_flag1 = np.array(list())
temp = lst1.value
for i in range(len(temp)):
	flag = True
	for ii in lst2.value:
		print(abs(temp[i] - ii))
		if( abs(temp[i] - ii) < 1.99):
			flag = False
			break
	if(flag):
		time_flag1 = np.append(time_flag1, time1[i] )
time_flag2 = np.array(list())
temp = lst2.value
print(temp)
for i in range(len(temp)):
	flag = True
	for ii in lst1.value:
		if( abs(temp[i] - ii) < 1.99):
			flag = False
			break
	if(flag):
		time_flag2 = np.append(time_flag2, time2[i] )
if(time_flag1.size != 0):
	print("to flag in obs1: " + str(time_flag1))
else:
	print("No time flags for obs1")
if(time_flag2.size != 0):
	print("to flag in obs2: " + str(time_flag2))
else:
	print("No time flags for obs2")



#Creation of new differenced ms

obs12 = root + "{0}-{1}-diff/{0}-{1}-diff.ms".format(obsid1, obsid2)
# Have to make a copy first
mset_diff = mset1.copy( obs12, deep = True)
# Open this copy
mset_diff = table( obs12, readonly = False)





#Get corrected data from which bad data is to be removed
data1 = mset1.getcol("CORRECTED_DATA")
data2 = mset2.getcol("CORRECTED_DATA")



#get indexes to be removed
#time
time_indexes1 = np.array(list())
times = mset1.getcol("TIME")
for i in time_flag1:
    time_indexes1 = np.append(time_indexes1,np.where(times == i))

time_indexes2 = np.array(list())
times = mset2.getcol("TIME")
for i in time_flag2:
    time_indexes2 = np.append(time_indexes2,np.where(times == i))

#antennas
ant_indexes1 = np.array(list())
ants = mset1.getcol("ANTENNA1")
for i in ants_to_flag:
    ant_indexes1 = np.append(ant_indexes1,np.where(ants == i))
ants = mset1.getcol("ANTENNA2")
for i in ants_to_flag:
    ant_indexes1 = np.append(ant_indexes1,np.where(ants == i))

ant_indexes2 = np.array(list())
ants = mset2.getcol("ANTENNA1")
for i in ants_to_flag:
    ant_indexes2 = np.append(ant_indexes2,np.where(ants == i))
ants = mset2.getcol("ANTENNA2")
for i in ants_to_flag:
    ant_indexes2 = np.append(ant_indexes2,np.where(ants == i))


print("time1:")
print(time_indexes1)

print("time2:")
print(time_indexes2)

#Removing from lists and ms
to_remove1 = np.unique(np.append(time_indexes1, ant_indexes1))
print(to_remove1.dtype)
to_remove1 = to_remove1.astype('int')
print(to_remove1.dtype)
#print(to_remove1)
to_remove2 = np.unique(np.append(time_indexes2, ant_indexes2))
print(to_remove2.dtype)
to_remove2 = to_remove2.astype('int')
print(to_remove2.dtype)
#NOTE: could remove data1 and just grab the information from the
# diff ms

data1 = np.delete(data1,to_remove1, 0)
data2 = np.delete(data2,to_remove2, 0)
mset_diff.removerows(to_remove1)
#NOTE: could remove data1 and just grab the information from the
# diff ms

#subtract the data
diff = data1 - data2

#put it into ms
mset_diff.putcol("CORRECTED_DATA", diff)




#match flags from ms2 to diff ms
flags = np.delete(mset1.getcol('FLAG'),to_remove1,0)
flags2 = np.delete(mset2.getcol('FLAG'),to_remove2,0)
num_rows=len(flags)
num_chans=len(flags[0])
asd = 0
#print( (flags == flags2).all() ) #tests for array equality
print("size of flag arrays:")
print("flag1: " + str(flags.shape))
print("flag2: " + str(flags2.shape))
for i in range(num_rows):
    for ii in range(num_chans):
        #if( flags2[i][ii][0] or flags[i][ii][0]):
        #if( flags2[i][ii][0] ):
        if( flags2[i][ii][0] and not flags[i][ii][0]):
            asd = asd + 1
            flags[i][ii][0] = True
            flags[i][ii][1] = True
            flags[i][ii][2] = True
            flags[i][ii][3] = True
print(str(asd) + " flags updated")

#put into diff ms
mset_diff.putcol("FLAG", flags)




flags = mset_diff.getcol("FLAG")


#for i in range(len(data1)):
#    for ii in range(len(data1[0])):
#        if (abs(diff[i][ii][0]) > 10**33):
#            print(diff[i][ii][0])
#            print(data1[i][ii][0])
#            print(data2[i][ii][0])
#            print(flags[i][ii][0])



mset1.close()
mset2.close()
mset_diff.close()

print("Subtraction Complete")
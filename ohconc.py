#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 19:23:09 2023

@author: maryamfatima
"""
import iris
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import h5py

OH2008 = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/OH/OH2008.nc')[0]
OH2009 = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/OH/OH2009.nc')[0]
OH2010 = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/OH/OH2010.nc')[0]
OH2011 = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/OH/OH2011.nc')[0]
OH2012 = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/OH/OH2012.nc')[0]

data08= OH2008.extract(iris.Constraint(time = lambda cell: cell.point.month == 3)) 
data09= OH2009.extract(iris.Constraint(time = lambda cell: cell.point.month == 11)) 
data10= OH2010.extract(iris.Constraint(time = lambda cell: cell.point.month == 3)) 
data11= OH2011.extract(iris.Constraint(time = lambda cell: cell.point.month == 3)) 
data12= OH2012.extract(iris.Constraint(time = lambda cell: cell.point.month == 3)) 

data08.units=("molecules cm-3")
level_heights = data08.coord('atmosphere_hybrid_height_coordinate').points


#%%
i=5
height=level_heights[i]
plt.style.use('classic')
data0=data10.extract(iris.Constraint(model_level_number = lambda cell: cell == i))
#data1= data0.collapsed('time', iris.analysis.MEAN)
conc=data0.data.data
conc=conc*1e-6
latitude = data0.coord('latitude').points
longitude=data0.coord('longitude').points
plt.figure(figsize=(21,16))
plt.rcParams.update({'font.size': 35})
plt.pcolor(longitude,latitude,conc)
#plt.clim(0,10)
cbar = plt.colorbar()
plt.xlim(0,360)
plt.ylim(-90,90)
cbar.set_label('OH Concentration (molecules cm$^{-3}$)')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
plt.title("OH Concentration (molecules cm$^{-3}$)\nin November 2010, %i m" % height)

#%%
for i in range (1,85):
    height=level_heights[i]
    plt.style.use('classic')
    removal_data_jan_50f=data08.extract(iris.Constraint(model_level_number = lambda cell: cell == i))
    coords = ['time']
    #removal_data_jan_50 = removal_data_jan_50f.collapsed(coords, iris.analysis.MEAN)
    conc=removal_data_jan_50f.data.data
    latitude = removal_data_jan_50f.coord('latitude').points
    longitude=removal_data_jan_50f.coord('longitude').points
    plt.figure(figsize=(21,16))
    plt.rcParams.update({'font.size': 35})
    plt.pcolor(longitude,latitude,conc)
    cbar = plt.colorbar()
    plt.xlim(0,360)
    plt.ylim(-90,90)
    plt.clim(0,5e6)
    cbar.set_label('$OH concentration in molecules per cm^3)')
    plt.title("OH concentration\n in March 2008, at %i m" % height)
    #plt.legend(loc='lower left',fontsize=8)
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    
    
    
    
    
    
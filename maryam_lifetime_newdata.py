#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 09:39:56 2023

@author: maryamfatima
"""

"""
lifetime, trying to do this quickly 
"""

import iris
import numpy as np
import matplotlib.pyplot as plt

CO14_data = iris.load_cube('/Volumes/Maryam backup/bsc_project/raw_data/CO14/CO142008.nc') #units mol/cm^3
CO14_prod = iris.load('/Volumes/Maryam backup/bsc_project/raw_data/production/CO14Production2010.nc')[0] #units kg/m2/s
CO14_loss = iris.load_cube('/Volumes/Maryam backup/bsc_project/raw_data/losses/2008Loss.nc')  #mol/gridcell/s
tropo = iris.load_cube('/Volumes/Maryam backup/bsc_project/raw_data/tropopause/2010Tropopause.nc')  #m
CO14_prod.aux_coords[0].rename('level_height')


CO14_data1=iris.load_cube('/Volumes/Maryam backup/bsc_project/raw_data/2010CO14.nc')
"""
uneven_levels = [-96.0, -72, -48, -24, -12, -6, 0, 6, 12, 24, 48, 72, 96]
vmin,vmax = -48,48
cs=plt.pcolormesh(lon,lat, data,cmap=cmap, transform=crs_new,vmin=vmin,vmax=vmax)
cbar=plt.colorbar(cs,boundaries= uneven_levels)
"""

#%%
months=['January','February','March','April','May','June','July','August','September','October','November','December']

#height at each grid cell 

h = CO14_data.coord('level_height').bounds[:,1]

###### Determining conversion factors #########
# mol to mass (kg) 
kg_mol = 1/30e-3 
#altitude of each cell
level_thick = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]

#CONVERSION FACTOR FOR PRODUCTION
prod_con = np.zeros([85,144,192])
prod_conv1 = kg_mol/(level_thick*10e6) #dividing by cell height and converting from m^3 to cm^3

for i in range(0,85):
    prod_con[i,:] = prod_conv1[i]

prod_conv = prod_con

#%%
#CONVERSION FACTOR FOR LOSS 

#finding box area
weights=iris.analysis.cartography.area_weights(CO14_data1)
area = weights[0,:,:,:]


#creating (85,144) array for height 
heights = np.zeros([85,144,192])
for i in range(0,85):
   heights[i,:] = level_thick[i]

volume = area * heights  
loss_conv = 1/(volume) 
print(volume)

#%%

conc0=CO14_data*1e6
loss0=CO14_loss*6.02e23*loss_conv

#%%
latitude = CO14_data.coord('latitude').points
altitude = CO14_data.coord('level_height').points
pressure = 1000 * np.exp( -altitude / 8.5e3)  
    
#%%
latitude = CO14_data.coord('latitude').points
longitude = CO14_data.coord('longitude').points
altitude = CO14_data.coord('level_height').points
pressure = 1000 * np.exp( -altitude / 8.5e3)

for j in range(1,13):
    index=j-1
    month=months[index]
    conc1= conc0.extract(iris.Constraint(time = lambda cell: cell.point.month == j))
    conc2= conc1.collapsed(['time'], iris.analysis.MEAN)
    conc3= conc2.collapsed(['longitude'], iris.analysis.SUM)
    loss1= loss0.extract(iris.Constraint(time = lambda cell: cell.point.month == j))
    loss2= loss1.collapsed(['time'], iris.analysis.MEAN)
    loss3= loss2.collapsed(['longitude'], iris.analysis.SUM)
    plt.style.use("classic")
    
    height=altitude[j]
    
    lifetime0=conc3/loss3
    lifetime1=lifetime0.data
    lifetime_day=lifetime1/60/60/24
    lifetime_month=lifetime_day/30
    plt.figure()
    plt.pcolor(longitude,latitude,lifetime_month)
    plt.gca().invert_yaxis()
    cbar = plt.colorbar()
    cbar.set_label('Lifetime of CO14 (months)')
    plt.title("Lifetime of CO14 (months)\n in %s 2008" % month)
    plt.xlim(0,360)
    plt.clim(0,12)
    plt.ylim(-90,90)
    #plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hPa)')
    plt.savefig('/Volumes/Maryam backup/bsc_project/lifetime_2008month{}.png'.format(j))   
    
#%%    
    
for i in range(1,13):
    index=i-1
    month=months[index]
    conc1= conc0.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    conc3= conc1.collapsed(['time','longitude'], iris.analysis.MEAN)
    loss1= loss0.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    loss3= loss1.collapsed(['time','longitude'], iris.analysis.MEAN)
    lifetime0=conc3/loss3
    lifetime1=lifetime0.data
    lifetime_day=lifetime1/60/60/24
    lifetime_month=lifetime_day/30
    plt.figure()
    plt.pcolor(latitude,pressure,lifetime_month)
    plt.gca().invert_yaxis()
    cbar = plt.colorbar()
    cbar.set_label('Lifetime of CO14 (months)')
    plt.title("Lifetime of CO14 (months)\n in %s 2008" % month)
    plt.xlim(-90,90)
    plt.clim(0,12)
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hPa)')
    plt.savefig('/Volumes/Maryam backup/bsc_project/lifetime_2010/lifetime2008month{}.png'.format(i))
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 22:40:56 2023

@author: ioanabalabasciuc
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import cftime
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs

#loading data

CO14_data = iris.load_cube('2010CO14.nc')

if not CO14_data.coord('latitude').has_bounds():
    CO14_data.coord('latitude').guess_bounds()
    CO14_data.coord('longitude').guess_bounds()
    

#%% Surface 

#CO14_data_surface = CO14_data.extract(iris.Constraint(model_level_number = lambda cell: cell == 1))

# coords_lt = ['longitude','time']
# data_mean = CO14_data_surface.collapsed(coords_lt,iris.analysis.MEAN, 
#                             weights = iris.analysis.cartography.area_weights(CO14_data_surface))
# plt.plot(data_mean.data.data)
# plt.ylabel('$^{14}$CO Concentration (molecule cm$^{-3}$)')
# plt.xlabel('Latitude box #')
# plt.title('Average Annual change with latitude for mean longitude')
# plt.show()


def concen_month(data,monthspec):
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN, 
                               weights = iris.analysis.cartography.area_weights(data_month))

    
    return data_mean


#creating 3d maps 
for m in range(1,13): 
    data_months = concen_month(CO14_data,m)
    plt.imshow(data_months.data.data.data)
    plt.title('14C concentration in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Altitude')
    plt.colorbar()
    plt.savefig('14C_con_long{}.png'.format(m))
    plt.show()
#in the console. This should be a 360 x 85 x 144 x 192. Because of the way these netCDF files



#%%
for j in range(1,10):
    #print('j',j)
    
    CO14_data_surface = CO14_data.extract(iris.Constraint(model_level_number = lambda cell: cell == j))
   
    surface = np.zeros(shape=(144,13))

    for m in range(1,13): 
        #print('m',m)
        data_months = concen_month(CO14_data_surface,m)
        surface[:,m] = data_months.data.data

    plt.plot(surface[:,1])
    plt.plot(surface[:,2])
    plt.plot(surface[:,3])
    plt.plot(surface[:,4])                              
    plt.plot(surface[:,5])
    plt.plot(surface[:,6])
    plt.plot(surface[:,7])
    plt.plot(surface[:,8])
    plt.plot(surface[:,9])
    plt.plot(surface[:,10])
    plt.plot(surface[:,11])
    plt.plot(surface[:,12])

    plt.ylabel('$^{14}$CO Concentration (molecule cm$^{-3}$)')
    plt.xlabel('Longitude box #')
    plt.title('Change with latitude for mean longitude, Cell= {}'.format(j))
    
    plt.legend(['1','2','3','4','5','6','7','8','9','10','11','12'],loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

    #total[:,j:j+12]
    
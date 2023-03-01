#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 17:16:32 2023

@author: ioanabalabasciuc
"""
import iris
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.colors as mcolors


CO14_08 = iris.load_cube('Concentration/CO142008.nc') #units mol/cm^3
CO14_09 = iris.load_cube('Concentration/CO142009.nc') #units mol/cm^3
CO14_10 = iris.load_cube('Concentration/CO142010.nc') #units mol/cm^3
CO14_11 = iris.load_cube('Concentration/CO142011.nc') #units mol/cm^3
CO14_12 = iris.load_cube('Concentration/CO142012.nc') #units mol/cm^3


tropo08 = iris.load_cube('Tropopause/2008Tropopause.nc')
tropo09 = iris.load_cube('Tropopause/2009Tropopause.nc')
tropo10 = iris.load_cube('Tropopause/2010Tropopause.nc')
tropo11 = iris.load_cube('Tropopause/2011Tropopause.nc')
tropo12 = iris.load_cube('Tropopause/2012Tropopause.nc')#m

#loading production files for all years 
CO14_prod08 = iris.load('Production/CO14Production2008.nc')[0] #units kg/m2/s
CO14_prod09 = iris.load('Production/CO14Production2009.nc')[0] #units kg/m2/s
CO14_prod10 = iris.load('Production/CO14Production2010.nc')[0] #units kg/m2/s
CO14_prod11 = iris.load('Production/CO14Production2011.nc')[0] #units kg/m2/s
CO14_prod12 = iris.load('Production/CO14Production2012.nc')[0] #units kg/m2/s

#loading loss files for all years 
CO14_loss08 = iris.load_cube('Loss/2008Loss.nc')  #mol/gridcell/s
CO14_loss09 = iris.load_cube('Loss/2009Loss.nc')  #mol/gridcell/s
CO14_loss10 = iris.load_cube('Loss/2010Loss.nc')  #mol/gridcell/s
CO14_loss11 = iris.load_cube('Loss/2011Loss.nc')  #mol/gridcell/s
CO14_loss12 = iris.load_cube('Loss/2012Loss.nc')  #mol/gridcell/s

tropo08 = iris.load_cube('Tropopause/2008Tropopause.nc')
tropo09 = iris.load_cube('Tropopause/2009Tropopause.nc')
tropo10 = iris.load_cube('Tropopause/2010Tropopause.nc')
tropo11 = iris.load_cube('Tropopause/2011Tropopause.nc')
tropo12 = iris.load_cube('Tropopause/2012Tropopause.nc')#m

vol_model08 = iris.load('Volume/2008Volume.nc')[0] #m3
vol_model09 = iris.load('Volume/2009Volume.nc')[0]#m3
vol_model10 = iris.load('Volume/2010Volume.nc')[0] #m3
vol_model11 = iris.load('Volume/2011Volume.nc')[0] #m3
vol_model12 = iris.load('Volume/2012Volume.nc')[0] #m3

#%%
kg_vol = []
for monthspec in range(1,13):
    data_month = CO14_10.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)
    col_month = vol_model08.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    
    CO14 = data_mean * col_month *1e6
    netsum = CO14.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).data 
    net_sum_new =  netsum*(1/6.02e23)*0.03
    kg_vol.append(net_sum_new)
    

print('Total_concentration', sum(kg_vol)/12)

#%% Converting Production from kg/m2/s to mol/m3/s
#height at each grid cell 
h = CO14_10.coord('level_height').bounds[:,1]

###### Determining conversion factors #########
# mol to mass (kg) 
kg_mol = 1/30e-3 
#altitude of each cell
level_thick = CO14_10.coord('level_height').bounds[:,1] -  CO14_10.coord('level_height').bounds[:,0]
#CONVERSION FACTOR FOR PRODUCTION
prod_con = np.zeros([85,144])
prod_conv1 = kg_mol/(level_thick*1e6) #dividing by cell height and converting from m^3 to cm^3

for i in range(0,85):
    prod_con[i,:] = prod_conv1[i]

prod_conv = prod_con

def concen_month(data,monthspec):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)

    return data_mean 

#%% Calcualting chemical lifetime  at  surface IN 2010 
for i in [1,7]:
    
    #obtaining loss for each month at surface 
    loss_month = CO14_loss11.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    loss_surface = loss_month.extract(iris.Constraint(model_level_number= lambda cell: cell == 1))
    rloss = loss_surface.collapsed('time',iris.analysis.MEAN)

    #obtaining cocentration for each month at surface 
    concen_month = CO14_11.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    concen_surface = concen_month.extract(iris.Constraint(model_level_number= lambda cell: cell == 1))
    rconc = concen_surface.collapsed('time',iris.analysis.MEAN).data
    
    rloss_array = np.array(rloss.data,dtype='float64')
    #converting production - loss into molecules/cm3/month
    loss_months1 = (rloss_array*6.02e23*(24*60*60))/(vol_model11.data[i-1,1,:,:]*1e6)
    
    #Calculating the chemical lifetime 
    data_months = rconc/loss_months1
    
    #extracting latitude at mindpoint of each cell in m 
    latitude = rloss.coord('latitude').points
    longitude = rloss.coord('longitude').points
    
    #creating uneven levels 
    uneven_levels = [10, 15,30, 60, 360]#720,360*3,4*360,5*360]
    cmap_rb = plt.get_cmap('Spectral')
    colors = cmap_rb(np.linspace(0, 1, len(uneven_levels) - 1))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    plt.figure()
    cs= plt.pcolormesh(longitude,latitude,data_months,cmap=cmap,norm=norm) # plotting latitude

    cbar= plt.colorbar(cs,ticks=uneven_levels,extend='max')
   
    cbar.set_label('Lifetime (day)')
    plt.title('Chemical lifetime of 14CO in 2011 #{} at 990mb'.format(i))
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')
    plt.savefig('ChemLoss/Surface/14C_chemloss11_{}.png'.format(i))
    plt.show()

#%% Calcualting chemical lifetime  at  550mb IN 2010 
for i in [1,7]:
    
    #obtaining loss for each month at surface 
    loss_month = CO14_loss11.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    loss_surface = loss_month.extract(iris.Constraint(model_level_number= lambda cell: cell == 29))
    rloss = loss_surface.collapsed('time',iris.analysis.MEAN)

    #obtaining cocentration for each month at surface 
    concen_month = CO14_11.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    concen_surface = concen_month.extract(iris.Constraint(model_level_number= lambda cell: cell == 29))
    rconc = concen_surface.collapsed('time',iris.analysis.MEAN).data
    
    rloss_array = np.array(rloss.data,dtype='float64')
    #converting production - loss into molecules/cm3/month
    loss_months1 = (rloss_array*6.02e23*(24*60*60))/(vol_model11.data[i-1,29,:,:]*1e6)
    
    #Calculating the chemical lifetime 
    data_months = rconc/loss_months1
    
    #extracting latitude at mindpoint of each cell in m 
    latitude = rloss.coord('latitude').points
    longitude = rloss.coord('longitude').points
    
    #creating uneven levels 
    uneven_levels = [10, 15,30, 60, 360]#720,360*3,4*360,5*360]
    cmap_rb = plt.get_cmap('Spectral')
    colors = cmap_rb(np.linspace(0, 1, len(uneven_levels) - 1))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    plt.figure()
    cs= plt.pcolormesh(longitude,latitude,data_months,cmap=cmap,norm=norm) # plotting latitude

    cbar= plt.colorbar(cs,ticks=uneven_levels,extend='max')
   
    cbar.set_label('Lifetime (day)')
    plt.title('Chemical lifetime of 14CO in 2011 #{} at 550mb'.format(i))
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')
    plt.savefig('ChemLoss/Surface/14C_chemloss550_11_{}.png'.format(i))
    plt.show()
    
    
#%% Calcualting chemical lifetime IN 2010 
for i in range(1,13):
    
    trop_month = tropo12.extract(iris.Constraint(time = lambda cell: cell.point.month == i))
    coords_lt = ['longitude','time']
    trop_mean = trop_month.collapsed(coords_lt,iris.analysis.MEAN)
   
    rconc = concen_month(CO14_12,i).data
    rloss = concen_month(CO14_loss12,i)  
    
    loss_months1 = (rloss.data*6.02e23*(24*60*60))/(vol_model12.data[i-1,:,:,0]*1e6)
    
    #calculating the lifetime 
    data_months = rconc/loss_months1

    #extracting latitude at mindpoint of each cell 
    latitude = rloss.coord('latitude').points
    altitude = (rloss.coord('level_height').points)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    trop_level = 1000 * np.exp( -trop_mean.data/ 8.5e3) # converting tropopause from m to hBar
    
    uneven_levels = [10, 15,30, 60, 360,3.6e3,3.6e5]
    cmap_rb = plt.get_cmap('Spectral')
    colors = cmap_rb(np.linspace(0, 1, len(uneven_levels) - 1))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, colors)

    plt.figure()
    cs= plt.pcolormesh(latitude,pressure,data_months,cmap=cmap,norm=norm) # plotting latitude
    cbar= plt.colorbar(cs,ticks=uneven_levels,extend='max')
    plt.plot(latitude,trop_level, color='black', label='Tropopause') # plot ropopause 
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 

    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar.set_label('Lifetime (day)')
    plt.title('Chemical lifetime in 2012 # %.0f ' %(i))#,#np.max(data_months)/360))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    plt.legend()
    plt.savefig('ChemLoss/ChemLoss2012/14C_chemloss{}.png'.format(i))
    plt.show()

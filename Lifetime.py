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
#loading data
CO14_data = iris.load_cube('2010CO14.nc') #units mol/cm^3

#loading production files for all years 
CO14_prod10 = iris.load('CO14Production2010.nc')[0] #units kg/m2/s
CO14_prod11 = iris.load('CO14Production2011.nc')[0] #units kg/m2/s
CO14_prod12 = iris.load('CO14Production2012.nc')[0] #units kg/m2/s

#loading loss files for all years 
CO14_loss10 = iris.load_cube('OHCO142010reaction.nc')  #mol/gridcell/s
CO14_loss11 = iris.load_cube('OHCO142011reaction.nc')  #mol/gridcell/s
CO14_loss12 = iris.load_cube('OHCO142012reaction.nc')  #mol/gridcell/s


tropo = iris.load_cube('2010Tropopause.nc')  #m

#%%
#height at each grid cell 
h = CO14_data.coord('level_height').bounds[:,1]

###### Determining conversion factors #########
# mol to mass (kg) 
kg_mol = 1/30e-3 
#altitude of each cell
level_thick = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]
print(np.shape(level_thick))
#CONVERSION FACTOR FOR PRODUCTION
prod_con = np.zeros([85,144])
prod_conv1 = kg_mol/(level_thick*10e6) #dividing by cell height and converting from m^3 to cm^3

for i in range(0,85):
    prod_con[i,:] = prod_conv1[i]

prod_conv = prod_con

#CONVERSION FACTOR FOR LOSS 

#finding box area
weights=iris.analysis.cartography.area_weights(CO14_loss10)
area = weights[0,:,:,0]
print(np.shape(area))
#creating (85,144) array for height 

heights = np.zeros([85,144])
for i in range(0,85):
   heights[i,:] = level_thick[i]
   


volume = area * heights  
print(volume[0,0])
print(volume[0,-1])
np.savetxt('volume.csv',volume*10e6)
loss_conv = 1/(volume*10e6) 

def concen_month(data,monthspec):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)

    return data_mean 

plt.plot(volume[0,:])
plt.yscale('log')

#%% Calcualting chemical lifetime 
for i in range(1,2):

    rloss = concen_month(CO14_loss12,i)
    rconc = concen_month(CO14_data,i)

    rloss.units = 'kg m-2 s-1'
    #converting production - loss into mole/cm3/month
    loss_months1 = (rloss.data*loss_conv)*6.02e23*(365*24*60*60)
    np.savetxt('lossJan.csv',rloss)
    np.savetxt('lossJanconv.csv', loss_months1)
    #obtaining concentration at month 
    data_months = np.divide(rconc.data,loss_months1)
    #extracting latitude at mindpoint of each cell 
    latitude = rloss.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 

    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    altitude = (rloss.coord('level_height').points)
    
    plt.figure(figsize = (25,12.5))
    plt.plot(latitude,loss_months1[65,:], label = 'Loss molecules/cm3/year')
    plt.plot(latitude,rconc.data[65,:],label='Concentration molecules/cm3')
    plt.legend()
    plt.ylim(bottom = 0.0)
    plt.show()
    
    # #defining the tropopause level - only works for pressure plot 
    # Tropopause = concen_month(tropo,i)    
    # trop_level = 1000 * np.exp( -Tropopause.data/ 8.5e3) # converting tropopause from m to hBar    
    
    # #using first order approximation to find pressure 
    # pressure = 1000 * np.exp( -altitude / 8.5e3)
    # plt.pcolor(latitude,pressure, data_months,vmax=10) # plotting latitude
    
    # plt.plot(latitude,trop_level, color='black', label='Tropopause') # plot ropopause 
    # plt.gca().invert_yaxis() # invert axis as max pressure at surface 
        
    # plt.yticks([10,100,200,300,400,500,600,700,800,900,000])

    # cbar = plt.colorbar()
    # cbar.set_label('Lifetime (year)')
    # plt.title('Chemical lifetime of 14CO in years #{}'.format(i))
    # plt.xlabel('Latitude')
    
    # plt.ylabel('Pressure (hbar)')
    # #plt.ylabel('Altitude (km)')
    # plt.legend(loc='lower left')
    # plt.savefig('ChemLoss2012/14C_chemloss{}.png'.format(i))
    # plt.show()
#%%
CO14_loss_surface = CO14_loss12.extract(iris.Constraint(model_level_number = lambda cell: cell == 85 ))

CO14_loss_surface_30S_60S = CO14_loss_surface.extract(iris.Constraint(latitude = lambda cell: 0<cell<60))
CO14_loss_surface_30N_60N = CO14_loss_surface.extract(iris.Constraint(latitude = lambda cell: -60<cell<0))
coords_30_60 = ['latitude','longitude']

CO14_loss_surface_30N_60N_mean = CO14_loss_surface_30N_60N.collapsed(coords_30_60,
                                                                 iris.analysis.MEAN,
                                                                 weights = iris.analysis.cartography.area_weights(CO14_loss_surface_30N_60N))

CO14_loss_surface_30S_60S_mean = CO14_loss_surface_30S_60S.collapsed(coords_30_60,
                                                                 iris.analysis.MEAN,
                                                                 weights = iris.analysis.cartography.area_weights(CO14_loss_surface_30S_60S))

plt.figure(figsize = (25,12.5))
plt.plot(CO14_loss_surface_30N_60N_mean.data.data *6.02e23,label='0N-60N')
plt.plot(CO14_loss_surface_30S_60S_mean.data.data*6.02e23, label='60S-0N')

plt.ylabel('$^{14}$CO Loss (molecules/s')
plt.xlabel('Days in 2010')
plt.title('$^{14}$CO Loss at surface for everyday in 2010')
plt.legend()
plt.ylim(bottom = 0.0)
plt.show()




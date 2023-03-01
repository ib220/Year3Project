# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:48:03 2023

@author: lpb20
"""

#You can always go ahead and make changes to colorbars and other plotting variables
# to try to ake the figure look either prettier, or  more informative.

import iris
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#loading data
CO14_08 = iris.load_cube('Concentration/CO142008.nc') #units mol/cm^3
CO14_09 = iris.load_cube('Concentration/CO142009.nc') #units mol/cm^3
CO14_10 = iris.load_cube('Concentration/CO142010.nc') #units mol/cm^3
CO14_11 = iris.load_cube('Concentration/CO142011.nc') #units mol/cm^3
CO14_12 = iris.load_cube('Concentration/CO142012.nc') #units mol/cm^3

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

tropo08 = iris.load_cube('Tropopause/2008Tropopause.nc')
tropo09 = iris.load_cube('Tropopause/2009Tropopause.nc')
tropo10 = iris.load_cube('Tropopause/2010Tropopause.nc')
tropo11 = iris.load_cube('Tropopause/2011Tropopause.nc')
tropo12 = iris.load_cube('Tropopause/2012Tropopause.nc')#m


def concen_month(data,monthspec):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)

    return data_mean 


#%%looking at height variation 

#plotting heights at centre of each cell 
level_heights = CO14_loss10.coord('level_height').points
plt.figure()
plt.plot(level_heights, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Altitude (m)')
plt.show()

#finding bounds for each cell - thickness non unifrom greater at top of model 
level_thickness = CO14_loss10.coord('level_height').bounds[:,1] -  CO14_loss10.coord('level_height').bounds[:,0]
plt.figure()
plt.plot(level_thickness, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Model Level Thickness (m)')
plt.show()
        
#%% #plotting 14CO concentration for each month and saving data
for m in range(1,13): 
     data_months = concen_month(CO14_12,m)
          
     #defining the tropopause level - only works for pressure plot 
     Tropopause = concen_month(tropo12,m)
     
     #extracting latitude at mindpoint of each cell 
     latitude = data_months.coord('latitude').points
     altitude = (data_months.coord('level_height').points)

     
     #using first order approximation to find pressure 
     pressure = 1000 * np.exp( -altitude / 8.5e3)
     trop_level = 1000 * np.exp( -Tropopause.data/ 8.5e3) # converting tropopause from m to hBar
     plt.pcolor(latitude,pressure,data_months.data, vmax=65) #plotting as function of pressure
     plt.plot(latitude,trop_level, color='red', label='Tropopause') # plot ropopause 
     plt.gca().invert_yaxis() # invert axis as max pressure at surface 
     
     plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

     cbar = plt.colorbar()
     cbar.set_label('$^{14}$CO Concentration (molecules cm$^{-3}$ )')
     plt.title('14CO concentration in month #{} in 2012'.format(m))
     plt.xlabel('Latitude')
     plt.ylabel('Pressure (hbar)')
     #plt.ylabel('Altitude (km)')
     plt.legend(loc='lower left')
     plt.savefig('ConcPlot/ConcPlot2012/14CO_conc_lat{}.png'.format(m))
     plt.show()
 
#%%#plotting 14CO production for each month and saving data
for z in range(1,13): 
    data_months = concen_month(CO14_prod12,z)
    Tropopause = concen_month(tropo12,z)    
    
    # mol to mass (kg) 
    kg_mol = 1/30e-3 
    #altitude of each cell
    level_thick = data_months.coord('atmosphere_hybrid_height_coordinate').bounds[:,1] -  data_months.coord('atmosphere_hybrid_height_coordinate').bounds[:,0]

    #CONVERSION FACTOR FOR PRODUCTION
    prod_con = np.zeros([85,144])
    prod_conv1 = kg_mol/(level_thick*1e6) #dividing by cell height and converting from m^3 to cm^3

    for i in range(0,85):
        prod_con[i,:] = prod_conv1[i]

    prod_conv = prod_con
    
    #saving data as nc file

    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    altitude = (data_months.coord('atmosphere_hybrid_height_coordinate').points)
    
    #defining the tropopause level - only works for pressure plot 

    trop_level = 1000 * np.exp( -Tropopause.data/ 8.5e3) # converting tropopause from m to hBar

    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    #converting from kg/m2/s to molecules/cm3/month
    prod_new = (data_months.data*prod_conv)*6.02e23*(30*24*60*60) # multiplied by avogardos for molecules 
    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,prod_new,vmax=7) #plotting as function of pressure
    plt.plot(latitude,trop_level, color='red', label='Tropopause') # plot ropopause 

    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Production (molecules cm$^{-3}$ month$^{-1}$) ')
    plt.title('14CO Production in month #{} in 2012'.format(z))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    plt.legend(loc='lower left')
    plt.savefig('ProdPlot/ProdPlot2012/14CO_prod{}.png'.format(z))
    plt.show()


#%%#plotting 14CO loss for each month and saving data
for y in range(1,13):
    data_months = concen_month(CO14_loss11,y)
    Tropopause = concen_month(tropo11,y)    

    loss_months1 = (data_months.data*6.02e23*(30*24*60*60))/(vol_model11.data[y-1,:,:,0]*1e6)

    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    altitude = (data_months.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    trop_level = 1000 * np.exp( -Tropopause.data/ 8.5e3) # converting tropopause from m to hBar
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)

    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,loss_months1)#,vmax=15)#plotting as function of pressure
    plt.plot(latitude,trop_level, color='red', label='Tropopause') # plot ropopause 
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
        
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Loss (molecules cm$^{-3}$ month$^{-1}$)')
    plt.title('14CO Loss in month #{} in 2011'.format(y))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    plt.legend(loc='lower left')
    plt.savefig('LossPlot/LossPlot2011/14C_loss{}.png'.format(y))
    plt.show()
    

#%% Calcualting production-loss
for i in range(1,13):
    rprod = concen_month(CO14_prod12,i)
    rloss = concen_month(CO14_loss12,i)
    Tropopause = concen_month(tropo12,i)    
    vol = vol_model12.data[i-1,:,:,0]*1e6
    
    kg_mol = 1/30e-3 
    #altitude of each cell
    level_thick = rprod.coord('atmosphere_hybrid_height_coordinate').bounds[:,1] -  rprod.coord('atmosphere_hybrid_height_coordinate').bounds[:,0]

    #CONVERSION FACTOR FOR PRODUCTION
    prod_con = np.zeros([85,144])
    prod_conv1 = kg_mol/(level_thick*1e6) #dividing by cell height and converting from m^3 to cm^3

    for j in range(0,85):
        prod_con[j,:] = prod_conv1[j]

    prod_conv = prod_con

    loss_months1 = (rloss.data*6.02e23*(30*24*60*60))/(vol)
    prod_months1 =  rprod.data*prod_conv *6.02e23*(30*24*60*60)# mol to mass (kg) 
     
    data_months = prod_months1 - loss_months1
    #extracting latitude at mindpoint of each cell 
    latitude = rloss.coord('latitude').points
    altitude = (rloss.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    trop_level = 1000 * np.exp( -Tropopause.data/ 8.5e3) # converting tropopause from m to hBar    
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude
    plt.pcolor(latitude,pressure,data_months,cmap='seismic',vmax=10,vmin=-10) #plotting as function of pressure
    plt.plot(latitude,trop_level, color='black', label='Tropopause') # plot ropopause 
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
        
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('P-L (molecules cm$^{-3}$ month$^{-1}$)')
    plt.title('(P-L) in 14CO in month #{} in 2012'.format(i))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    plt.legend(loc='lower left')
    plt.savefig('NetPlot/NetPlot2012/14C_net_cont{}.png'.format(i))
    plt.show()


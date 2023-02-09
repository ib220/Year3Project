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
CO14_data = iris.load_cube('2010CO14.nc') #units mol/cm^3
CO14_prod = iris.load('CO14Production2010.nc')[0] #units kg/m2/s
CO14_loss = iris.load_cube('OHCO14reaction.nc')  #mol/gridcell/s

#%%
#unit used for production and loss will be mol/cm^3/s

#radius of the earth 
re= 6378.1e3
#height at each grid cell 
h = CO14_data.coord('level_height').bounds[:,1]

###### Determining conversion factors #########
# mol to mass (kg) 
kg_mol = 1/30e-3 
#altitude of each cell
level_thick = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]

#CONVERSION FACTOR FOR PRODUCTION
prod_con = np.zeros([85,144])
prod_conv1 = kg_mol/(level_thick*10e6) #dividing by cell height and converting from m^3 to cm^3

for i in range(0,85):
    prod_con[i,:] = prod_conv1[i]

prod_conv = prod_con

#CONVERSION FACTOR FOR LOSS 

#finding box area
weights=iris.analysis.cartography.area_weights(CO14_loss)
area = weights[0,:,:,0]


#creating (85,144) array for height 
heights = np.zeros([85,144])
for i in range(0,85):
   heights[i,:] = level_thick[i]

volume = area * heights  
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

#%%
#looking at height variation 

#plotting heights at centre of each cell 
level_heights = CO14_data.coord('level_height').points
plt.figure()
plt.plot(level_heights, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Altitude (m)')
plt.show()

#finding bounds for each cell - thickness non unifrom greater at top of model 
level_thickness = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]
plt.figure()
plt.plot(level_thickness, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Model Level Thickness (m)')
plt.show()
        
#%% #plotting 14CO concentration for each month and saving data

for m in range(1,13): 
     data_months = concen_month(CO14_data,m)
     #saving data as nc file
     iris.save(data_months,'ConcData/14COconc{}.nc'.format(m) )
     
     #extracting latitude at mindpoint of each cell 
     latitude = data_months.coord('latitude').points
     #extracting altitude at mindpoint of each cell  in m 
     
     ##### NOTE #######
     #height in concentration is level_height 
     #height in prod and loss is atmosphere_hybrid_height_coordinate
     #altitude = (data_months.coord('level_height').points)
     altitude = (data_months.coord('level_height').points)
     
     #defining the tropopause level - only works for pressure plot 
     Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
     
     #using first order approximation to find pressure 
     pressure = 1000 * np.exp( -altitude / 8.5e3)
     #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
     plt.pcolor(latitude,pressure,data_months.data, vmax=45) #plotting as function of pressure
     plt.gca().invert_yaxis() # invert axis as max pressure at surface 
     
     plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
     
     plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

     cbar = plt.colorbar()
     cbar.set_label('$^{14}$CO Concentration (moles cm$^{-3}$ )')
     plt.title('14CO concentration in month #{}'.format(m))
     plt.xlabel('Latitude')
     plt.ylabel('Pressure (hbar)')
     #plt.ylabel('Altitude (km)')
     plt.legend(loc='lower left')
     plt.savefig('ConcPlot/14C_conc_long{}.png'.format(m))
     plt.show()
 
#%%#plotting 14CO production for each month and saving data
for z in range(1,2): 
    data_months = concen_month(CO14_prod,m)
    #saving data as nc file
    iris.save(data_months,'ProdData/14COprod.nc')
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (data_months.coord('atmosphere_hybrid_height_coordinate').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    m = (data_months.data*prod_conv)*6.02e23*2.628e6 # multiplied by avogardos for molecules 
    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,m) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Production (molecules cm$^{-3}$ month$^{-1}$) ')
    #plt.title('14CO Production in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('ProdPlot/14C_prod.png')
    plt.show()


#%%#plotting 14CO loss for each month and saving data

for m in range(1,13): 
    data_months = concen_month(CO14_loss,m)
    #saving data as nc file
    iris.save(data_months,'LossData/14COloss{}.nc'.format(m) )
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (data_months.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude
    loss = data_months.data*loss_conv*6.02e23*2.628e6  # multiplied by avogardos for molecules 
    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,loss)# vmax=1)#plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Loss (molecules cm$^{-3}$ month$^{-1}$)')
    plt.title('14CO Loss in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('LossPlot/14C_loss_long{}.png'.format(m))
   
    plt.show()
    

#%% Looking at seasonality 
loss_jan = iris.load_cube('ConcData/14COconc1.nc') #units mol/cm^3
loss_july = iris.load_cube('ConcData/14COconc7.nc')

CO14seasnoaldiff = 100 * (loss_jan - loss_july) /  loss_jan 


latitude = CO14seasnoaldiff.coord('latitude').points
altiude = CO14seasnoaldiff.coord('level_height').points
pressure = 1000 * np.exp( -altiude / 8.5e3)
CO14seasnoaldiff_data = CO14seasnoaldiff.data

plt.figure()
plt.contourf(latitude,pressure,CO14seasnoaldiff_data, levels = np.linspace(-100,100,300), cmap = 'seismic')
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label('% change between January and July in $^{14}$CO Loss')
plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])
plt.xlabel('Latitude')
plt.ylabel('Pressure (hPa)')

#%% Calcualting production-loss
for i in range(1,13):
    rprod =  iris.load('ProdData/14COprod.nc')[0]
    rloss = iris.load_cube('LossData/14COloss{}.nc'.format(i))
    rloss.units = 'kg m-2 s-1'
    #print(  rloss)
    data_months = (((rprod.data*prod_conv  - rloss.data*loss_conv))*6.02e23)*2.628e6
    #extracting latitude at mindpoint of each cell 
    latitude = rloss.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (rloss.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    plt.pcolor(latitude,pressure,data_months)#, vmin=-1) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('P-L (molecules cm$^{-3}$ month$^{-1}$)')
    plt.title('(P-L)/P in 14CO in month #{}'.format(i))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('NetPlot/14C_net{}.png'.format(i))
    plt.show()
    





#CONVERSION FACTOR FOR PRODUCTION
prod_con = np.zeros([85,144])
prod_conv1 = kg_mol/(level_thick*10e6) #dividing by cell height and converting from m^3 to cm^3

for i in range(0,85):
    prod_con[i,:] = prod_conv1[i]

prod_conv = prod_con

#CONVERSION FACTOR FOR LOSS 

#finding box area
weights=iris.analysis.cartography.area_weights(CO14_loss)
area = weights[0,:,:,0]


#creating (85,144) array for height 
heights = np.zeros([85,144])
for i in range(0,85):
   heights[i,:] = level_thick[i]

volume = area * heights  
print(volume)
loss_conv = 1/(volume*10e6) 

def concen_month(data,monthspec, prod=0, loss=0):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)

    return data_mean 

#%%
#looking at height variation 

#plotting heights at centre of each cell 
level_heights = CO14_data.coord('level_height').points
plt.figure()
plt.plot(level_heights, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Altitude (m)')
plt.show()

#finding bounds for each cell - thickness non unifrom greater at top of model 
level_thickness = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]
plt.figure()
plt.plot(level_thickness, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Model Level Thickness (m)')
plt.show()
        
#%% #plotting 14CO concentration for each month and saving data

for m in range(1,13): 
     data_months = concen_month(CO14_data,m)
     #saving data as nc file
     iris.save(data_months,'ConcData/14COconc{}.nc'.format(m) )
     
     #extracting latitude at mindpoint of each cell 
     latitude = data_months.coord('latitude').points
     #extracting altitude at mindpoint of each cell  in m 
     
     ##### NOTE #######
     #height in concentration is level_height 
     #height in prod and loss is atmosphere_hybrid_height_coordinate
     #altitude = (data_months.coord('level_height').points)
     altitude = (data_months.coord('level_height').points)
     
     #defining the tropopause level - only works for pressure plot 
     Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
     
     #using first order approximation to find pressure 
     pressure = 1000 * np.exp( -altitude / 8.5e3)
     #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
     plt.pcolor(latitude,pressure,data_months.data, vmax=45) #plotting as function of pressure
     plt.gca().invert_yaxis() # invert axis as max pressure at surface 
     
     plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
     
     plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

     cbar = plt.colorbar()
     cbar.set_label('$^{14}$CO Concentration (moles cm$^{-3}$ )')
     plt.title('14CO concentration in month #{}'.format(m))
     plt.xlabel('Latitude')
     plt.ylabel('Pressure (hbar)')
     #plt.ylabel('Altitude (km)')
     plt.legend(loc='lower left')
     plt.savefig('ConcPlot/14C_conc_long{}.png'.format(m))
     plt.show()
 
#%%#plotting 14CO production for each month and saving data
for m in range(1,12): 
    data_months = concen_month(CO14_prod,m)
    #saving data as nc file
    iris.save(data_months,'ProdData/14COprod{}.nc'.format(m) )
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (data_months.coord('atmosphere_hybrid_height_coordinate').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    m = (data_months.data*prod_conv)*6.02e23*2.628e6 # multiplied by avogardos for molecules 
    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,m) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Production (mol cm$^{-3}$ s${-1})$ ')
    #plt.title('14CO Production in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('ProdPlot/14C_prod_long{}.png'.format(m))
    plt.show()


#%%#plotting 14CO loss for each month and saving data

for m in range(1,13): 
    data_months = concen_month(CO14_loss,m,loss=1)
    #saving data as nc file
    iris.save(data_months,'LossData/14COloss{}.nc'.format(m) )
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (data_months.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude
    loss = data_months.data*loss_conv*6.02e23*2.628e6  # multiplied by avogardos for molecules 
    #multipled by number of seconds in a month to get molecules/cm3/months
    plt.pcolor(latitude,pressure,loss)#, vmin=0, vmax=8e-31 ) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Loss (moles cm$^{-3}$ s$^{-1}$)')
    plt.title('14CO Loss in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('LossPlot/14C_loss_long{}.png'.format(m))
   
    plt.show()
    

#%% Looking at seasonality 
loss_jan = iris.load_cube('ConcData/14COconc1.nc') #units mol/cm^3
loss_july = iris.load_cube('ConcData/14COconc7.nc')

CO14seasnoaldiff = 100 * (loss_jan - loss_july) /  loss_jan 


latitude = CO14seasnoaldiff.coord('latitude').points
altiude = CO14seasnoaldiff.coord('level_height').points
pressure = 1000 * np.exp( -altiude / 8.5e3)
CO14seasnoaldiff_data = CO14seasnoaldiff.data

plt.figure()
plt.contourf(latitude,pressure,CO14seasnoaldiff_data, levels = np.linspace(-100,100,300), cmap = 'seismic')
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label('% change between January and July in $^{14}$CO Loss')
plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])
plt.xlabel('Latitude')
plt.ylabel('Pressure (hPa)')

#%% Calcualting P-L 
for i in range(1,13):
    rprod =  iris.load('ProdData/14COprod{}.nc'.format(i))[0]
    rloss = iris.load_cube('LossData/14COloss{}.nc'.format(i))
    rloss.units = 'kg m-2 s-1'
    #print(  rloss)
    data_months = (((rprod.data*prod_conv  - rloss.data*loss_conv))*6.02e23)*2.628e6
    #extracting latitude at mindpoint of each cell 
    latitude = rloss.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in m 
    
    ##### NOTE #######
    #height in concentration is level_height 
    #height in prod and loss is atmosphere_hybrid_height_coordinate
    #altitude = (data_months.coord('level_height').points)
    altitude = (rloss.coord('level_height').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    plt.pcolor(latitude,pressure,data_months) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('% (P-L)/P')
    plt.title('(P-L)/P in 14CO in month #{}'.format(i))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
    #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('NetPlot/14C_net{}.png'.format(i))
    plt.show()
    

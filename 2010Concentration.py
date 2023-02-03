#You can always go ahead and make changes to colorbars and other plotting variables
# to try to ake the figure look either prettier, or  more informative.

import iris
import numpy as np
import matplotlib.pyplot as plt

#loading data
CO14_data = iris.load_cube('2010CO14.nc')
CO14_prod = iris.load('CO14Production2010.nc')[0]
#CO14_loss = iris.load_cube('OHCO14reaction.nc')

#%%
#looking at height variation 

#plotting heights at centre of each cell 
level_heights = CO14_data.coord('level_height').points
plt.figure()
plt.plot(level_heights, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Altitude (m)')

#finding bounds for each cell - thickness non unifrom greater at top of model 
level_thickness = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]
plt.figure()
plt.plot(level_thickness, 'x')
plt.xlabel('Model Level Height')
plt.ylabel('Model Level Thickness (m)')


#%%


def concen_month(data,monthspec):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)
    
    return data_mean



#plotting 14CO concentration for each month and saving data
for m in range(1,13): 
    data_months = concen_month(CO14_data,m)
    #saving data as nc file
    iris.save(data_months,'ConcData/14COcon{}.nc'.format(m) )
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in km 
    altitude = (data_months.coord('level_height').points)/1000
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)/1000
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    plt.pcolor(latitude,pressure,data_months.data) #plotting presusre 
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause) # plot approx tropopause only for pressure
    
    #plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Concentration (molecules cm$^{-3}$)')
    plt.title('14C concentration in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Altitude (km)')
    plt.savefig('ConcPlot/14C_con_long{}.png'.format(m))
    plt.show()

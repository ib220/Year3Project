
#You can always go ahead and make changes to colorbars and other plotting variables
# to try to ake the figure look either prettier, or  more informative.

import iris
import numpy as np
import matplotlib.pyplot as plt

#loading data
CO14_data = iris.load_cube('2010CO14.nc') #units mol/cm^3
CO14_prod = iris.load('CO14Production2010.nc')[0] #units kg/m2/s
CO14_loss = iris.load_cube('OHCO14reaction.nc')  #mol/gridcell/s
#unit used for production and loss will be mol/cm^3/s

#Determining conversion factors
# mol to mass (kg) 
kg_mol = 1/30e-3 

#altitude of each cell
level_thick = CO14_data.coord('level_height').bounds[:,1] -  CO14_data.coord('level_height').bounds[:,0]
#longitude of each cell
long_len = CO14_data.coord('longitude').bounds[:,1] -  CO14_data.coord('longitude').bounds[:,0]
#latitude of each cell 
lat_len = CO14_data.coord('latitude').bounds[:,1] -  CO14_data.coord('latitude').bounds[:,0]
#longitute and latitude  box thickness is the same 

#conversion factor for production 
prod_conv = kg_mol/(level_thick * 10e-6) #dividing by cell height and converting from m^3 to cm^3
#conversion factor for loss 
loss_conv = 1/(lat_len[0]*long_len[0]*level_thick)

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

##### NOTE #######
#height in concentration is level_height 
#height in prod and loss is atmosphere_hybrid_height_coordinate

#%%


def concen_month(data,monthspec, prod=0, loss=0):
    """
    Function collpases cube over longitude and time, 
    returning mean 14CO concentration
    """
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)
    
    if prod==1:
        data_mean_new = prod_conv * data_mean
        return data_mean_new 
    else: 
        return data_mean 
    
    if loss==1: 
        data_mean_new = loss_conv * data_mean
        return data_mean_new 
    else: 
        return data_mean 
        
    

#plotting 14CO concentration for each month and saving data
for m in range(1,2): 
    data_months = concen_month(CO14_prod,m, prod=1)
    #saving data as nc file
    iris.save(data_months,'ProdData/14COprod{}.nc'.format(m) )
    
    #extracting latitude at mindpoint of each cell 
    latitude = data_months.coord('latitude').points
    #extracting altitude at mindpoint of each cell  in km 
   # altitude = (data_months.coord('level_height').points)
    altitude = (data_months.coord('atmosphere_hybrid_height_coordinate').points)
    
    #defining the tropopause level - only works for pressure plot 
    Tropopause = (300 - 215 * (np.cos((np.pi * latitude)/180)) ** 2)
    
    #using first order approximation to find pressure 
    pressure = 1000 * np.exp( -altitude / 8.5e3)
    
    #plt.pcolor(latitude,altitude,data_months.data) # plotting latitude 
    plt.pcolor(latitude,pressure,data_months.data) #plotting as function of pressure
    plt.gca().invert_yaxis() # invert axis as max pressure at surface 
    
    plt.plot(latitude,Tropopause, color='red', label='Tropopause') # plot approx tropopause only for pressure
    
    plt.yticks([10,100,200,300,400,500,600,700,800,900,1000])

    cbar = plt.colorbar()
    cbar.set_label('$^{14}$CO Production (molecules cm$^{-3}$)')
    plt.title('14CO production in month #{}'.format(m))
    plt.xlabel('Latitude')
    plt.ylabel('Pressure (hbar)')
   #plt.ylabel('Altitude (km)')
    plt.legend(loc='lower left')
    plt.savefig('ProdPlot/14C_prod_long{}.png'.format(m))
    plt.show()






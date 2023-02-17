#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:20:34 2023

@author: ioanabalabasciuc
"""
import iris
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#loading data
CO14_data = iris.load_cube('2010CO14.nc') #units mol/cm^3
CO14_prod = iris.load('CO14Production2010.nc')[0] #units kg/m2/s
CO14_loss = iris.load_cube('OHCO14reaction.nc')  #mol/gridcell/s
tropo = iris.load_cube('2010Tropopause.nc')  #m
CO14_prod.aux_coords[0].rename('level_height')

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
area_2d = weights[0,:,:,0]
area_3d = weights[0,:,:,:]
print(np.shape(area_3d))

#creating (85,144) array for height 
heights = np.zeros([85,144])
for i in range(0,85):
   heights[i,:] = level_thick[i]

volume = area_2d * heights  
loss_conv = 1/(volume*10e6) 
#%% FInding net loss and produciton 
CO14_prod_new = CO14_prod * iris.analysis.cartography.area_weights(CO14_prod)
netprod = CO14_prod_new.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
CO14_loss_new = CO14_loss * 30e-3
netloss = CO14_loss_new.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
print('net prod', netprod)
print('net loss', netloss)

#%%
#finding model level number above which data in tropopause 
def trop_barrier(tropo,data):
    for t in range(1,13):
        print('Month creating model level number', t)
        Tropopause = tropo.extract(iris.Constraint(time = lambda cell: cell.point.month == 1))   
        #finding if above tropopause 
        #defining array of latitudes for each tropopause value 
        lat_trop = Tropopause.coord('latitude').points
        long_trop = Tropopause.coord('longitude').points
        #model level number above tropopause for each latitude 
        heights = np.zeros([145,193])
        lat_numb=0
        for lat in lat_trop:
            lat_numb += 1
            long_numb=0
            for long in long_trop:
                long_numb +=1
                #obtaining tropopause height at latitude point 
                lat_constraint = iris.Constraint(latitude = lambda cell: cell==lat)

                long_constraint = iris.Constraint(longitude = lambda cell: cell==long)
            
                trop_val = Tropopause.extract(lat_constraint & long_constraint).data
            
                top_box = data.coord('level_height').bounds[:,1]

                #obtaining  model level number above tropopause for each latitude 
                box_number =0 
                for h in range(0,85): 
                    if top_box[h] < trop_val:
                        box_number = h+1
                        heights[int(lat_numb),int(long_numb)] = box_number

        np.savetxt('Barrier/barrier{}.csv'.format(t),heights , delimiter=',')
    
trop_barrier(tropo,CO14_data) 
#%% Creating list with mean values for each quadrant 
mean_prod_TS = []
mean_prod_TN = [] 
mean_prod_SS = [] 
mean_prod_SN = [] 

mean_loss_TS = []
mean_loss_TN = [] 
mean_loss_SS = [] 
mean_loss_SN = [] 

#%% Creating monthly 4 quadrants mask 
#finding model level number at which data in tropopause for each month 
for mon in range(1,13):
    print('Month # {} contracts',mon)
    # obtaining production and loss data in Jan 
    Jan_prod = CO14_prod_new.extract(iris.Constraint(time = lambda cell: cell.point.month == mon)).data* 3.15e7
    #convert to kg/s for kg\ms\s
    Jan_loss = CO14_loss_new.extract(iris.Constraint(time = lambda cell: cell.point.month == mon)).data* 3.15e7
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/barrier{}.csv'.format(mon),delimiter=',', unpack=1)
    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        for longi in longitudes:
            angle = (180/144)*lati
            for h in range(0,85):
                if h >= barrier_height[longi+1,lati+1] and angle>90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi+1,lati+1] and angle<90:
                    month_topN[h,lati,longi]=1
                if h <= barrier_height[longi+1,lati+1] and angle<90:
                    month_botN[h,lati,longi]=1
                if h <= barrier_height[longi+1,lati+1] and angle>90:
                    month_botS[h,lati,longi]=1
    
    month_topS
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topS
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topS
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
#%% creatiing plots of averages for loss and production
month_numbers = np.arange(1,1)
plt.figure()
plt.title('CO14 Production')
plt.plot(month_numbers,mean_prod_TS, 'x',label='Trop South')
plt.plot(month_numbers,mean_prod_TN, 'x',label='Trop North')
plt.plot(month_numbers,mean_prod_SS, 'x',label='Strat South')
plt.plot(month_numbers,mean_prod_SN, 'x',label='Strat North')
#plt.xticks([']','F','M','A','M','J','J','A','S','O','N','D'])
plt.xlabel('Month')
plt.ylabel('CO14 Production (kg year${-1}$)')
plt.legend()
plt.savefig('Average_Prod.pdf')
plt.show()

plt.figure()
plt.title('CO14 Loss')
plt.plot(month_numbers,mean_loss_TS, 'x',label='Trop South')
plt.plot(month_numbers,mean_loss_TN, 'x',label='Trop North')
plt.plot(month_numbers,mean_loss_SS, 'x',label='Strat South')
plt.plot(month_numbers,mean_loss_SN, 'x',label='Strat North')
#plt.xticks(['J','F','M','A','M','J','J','A','S','O','N','D'])
plt.xlabel('Month')
plt.ylabel('CO14 Loss (kg year${-1}$)')
plt.legend()
plt.savefig('Average_Loss.pdf')
plt.show()

print(sum(mean_prod_TS)+sum(mean_prod_TN) + sum(mean_prod_SS) + sum(mean_prod_SS) + sum(mean_prod_SN) )
print(sum(mean_loss_TS)+sum(mean_loss_TN) + sum(mean_loss_SS) + sum(mean_prod_SS) + sum(mean_loss_SN) )
print(sum(mean_prod_TS))


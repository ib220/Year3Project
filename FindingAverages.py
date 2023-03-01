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
vol1_model2 = iris.load('Volume/2012Volume.nc')[0] #m3

#%% Finding converstion factor for production from kg/m2/s to mol/cm^3/s

#height of each grid cell 
h = CO14_loss10.coord('level_height').bounds[:,1]

# mol to mass (kg) 
kg_mol = 1/30e-3 
#altitude of each cell
level_thick = CO14_loss10.coord('level_height').bounds[:,1] -  CO14_loss10.coord('level_height').bounds[:,0]

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
    data_month = data.extract(iris.Constraint(time = lambda cell: cell.point.month == monthspec))
    coords_lt = ['longitude','time']
    data_mean = data_month.collapsed(coords_lt,iris.analysis.MEAN)

    return data_mean 

#%%Finding net loss and production for the all the years 

#for 2010 
#finding production in kg/s
CO14_prod_rate10 = CO14_prod10 * iris.analysis.cartography.area_weights(CO14_prod10)
#finding total production rate in kg/s for a year
netprod10 = CO14_prod_rate10.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7

#converting loss from mol/s to kg/s
CO14_loss_kg10 = CO14_loss10 * 30e-3
#finding total production rate in kg/year 
netloss10 = CO14_loss_kg10.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7

#repeating for rest of years
 
#for 2008 
CO14_prod_rate08 = CO14_prod08 * iris.analysis.cartography.area_weights(CO14_prod08)
netprod08 = CO14_prod_rate08.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
CO14_loss_kg08 = CO14_loss08 * 30e-3
netloss08 = CO14_loss_kg08.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
                                      
#for 2009
CO14_prod_rate09 = CO14_prod09 * iris.analysis.cartography.area_weights(CO14_prod09)
netprod09 = CO14_prod_rate09.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
CO14_loss_kg09 = CO14_loss09 * 30e-3
netloss09 = CO14_loss_kg09.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7

#for 2011
CO14_prod_rate11 = CO14_prod11 * iris.analysis.cartography.area_weights(CO14_prod11)
netprod11 = CO14_prod_rate11.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
CO14_loss_kg11 = CO14_loss11 * 30e-3
netloss11 = CO14_loss_kg11.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7

#for 2012
CO14_prod_rate12 = CO14_prod12 * iris.analysis.cartography.area_weights(CO14_prod12)
netprod12 = CO14_prod_rate12.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
CO14_loss_kg12 = CO14_loss12 * 30e-3
netloss12 = CO14_loss_kg12.collapsed(['latitude','longitude','model_level_number'],iris.analysis.SUM).collapsed('time',iris.analysis.MEAN).data * 3.15e7
#%%

plt.plot([2008,2009,2010,2011,2012],[netprod08,netprod09,netprod10,netprod11,netprod12],'x',label='Loss',color='red')
plt.plot([2008,2009,2010,2011,2012],[netloss08,netloss09,netloss10,netloss11,netloss12],'x',label='Production',color='blue')

plt.xticks([2008,2009,2010,2011,2012])
plt.xlabel('Year')
plt.ylabel('Loss or Production (kg year ${-1}$)')
plt.title('Annual Loss and Production')
plt.grid()
plt.legend()
plt.show()

plt.plot([2008,2009,2010,2011,2012],[netprod08- netloss08,netprod09 -netloss09,netprod10-netloss10,netprod11-netloss11,netprod12-netloss12],'x',label='P-L',color='purple')
plt.xticks([2008,2009,2010,2011,2012])
plt.xlabel('Year')
plt.ylabel('Production  - Loss (kg year ${-1}$)')
plt.title('Annual Production - Loss ')
plt.grid()
plt.show()


#%%
#finding model level number above which data in tropopause 
def trop_barrier(tropo,data):
    for t in range(0,12):
        print('Month creating model level number', t)
        Tropopause = tropo.extract(iris.Constraint(time = lambda cell: cell.point.month == 1))   
        #finding if above tropopause 
        #defining array of latitudes for each tropopause value 
        lat_trop = np.arange(0,144)
        long_trop= np.arange(0,192)
        #model level number above tropopause for each latitude 
        heights = np.zeros([144,192])
       
        for lat in lat_trop:
            for long in long_trop:
                lat_trop_val = Tropopause.coord('latitude').points[int(lat)]
                long_trop_val = Tropopause.coord('longitude').points[int(long)]
                #obtaining tropopause height at latitude point 
                lat_constraint = iris.Constraint(latitude = lambda cell: cell==lat_trop_val)

                long_constraint = iris.Constraint(longitude = lambda cell: cell== long_trop_val)
            
                trop_val = Tropopause.extract(lat_constraint & long_constraint).data
            
                top_box = data.coord('atmosphere_hybrid_height_coordinate').bounds[:,1]

                #obtaining  model level number above tropopause for each latitude 
                box_number =0 
                height = np.arange(0,84)
                for h in height: 
                    if top_box[h] < trop_val:
                        box_number = h+1
                        heights[int(lat),int(long)] = box_number

        np.savetxt('Barrier/Barrier2009/barrier{}.csv'.format(t),heights , delimiter=',')
    
trop_barrier(tropo09,CO14_prod09) 
#%% Creating list with mean values for each quadrant 
mean_prod_TS = []
mean_prod_TN = [] 
mean_prod_SS = [] 
mean_prod_SN = [] 

mean_loss_TS = []
mean_loss_TN = [] 
mean_loss_SS = [] 
mean_loss_SN = [] 

mean_pl_TS = []
mean_pl_TN = [] 
mean_pl_SS = [] 
mean_pl_SN = [] 

#%% Creating monthly 4 quadrants mask 2008
#finding model level number at which data in tropopause for each month 
print('2008')
for mon in range(1,13):
    print('Month # {} contracts'.format(mon))
    # obtaining production and loss data in Jan 
    Jan_prod2 = CO14_prod08.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_prod1 = Jan_prod2.data * 3.15e7
    Jan_prod = Jan_prod1  * iris.analysis.cartography.area_weights(CO14_prod10)[mon-1,:,:,:]
    #convert to kg/s from kg\ms\s
    Jan_loss2 = CO14_loss08.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_loss1 = Jan_loss2.collapsed('time',iris.analysis.MEAN)
    Jan_loss = Jan_loss1.data * 3.15e7 * 30e-3
    #calculating net production and loss 
    Jan_pl = Jan_prod - Jan_loss 
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/Barrier2008/barrier{}.csv'.format(mon-1),delimiter=',', unpack=1)

    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        angle = (180/144)*lati
        for longi in longitudes:
            for h in range(0,85):
                if h >= barrier_height[longi,lati] and angle>=90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi,lati] and angle<90:
                    month_topN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle<90:
                    month_botN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle>=90:
                    month_botS[h,lati,longi]=1
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topN
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topN
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 
    
    Jan_pl_SS = Jan_pl * month_topS
    Jan_pl_SN = Jan_pl * month_topN
    Jan_pl_TS = Jan_pl * month_botS 
    Jan_pl_TN = Jan_pl * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)#/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)#/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)#/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)#/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)#/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)#/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)#/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)#)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
    #finding average over each quadrant for loss
    JanPLSS_avg = np.sum(Jan_pl_SS)#/np.sum(month_topS)
    mean_pl_SS.append(JanPLSS_avg)
    JanPLSN_avg = np.sum(Jan_pl_SN)#/np.sum(month_topN)
    mean_pl_SN.append(JanPLSN_avg)
    JanPLTN_avg = np.sum(Jan_pl_TN)#/np.sum(month_botN)
    mean_pl_TN.append(JanPLTN_avg)
    JanPLTS_avg = np.sum(Jan_pl_TS)#)/np.sum(month_botS)
    mean_pl_TS.append(JanPLTS_avg)
#%% Creating monthly 4 quadrants mask 2009
#finding model level number at which data in tropopause for each month 
print('2009')
for mon in range(1,13):
    print('Month # {} contracts'.format(mon))
    # obtaining production and loss data in Jan 
    Jan_prod2 = CO14_prod09.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_prod1 = Jan_prod2.data * 3.15e7
    Jan_prod = Jan_prod1  * iris.analysis.cartography.area_weights(CO14_prod10)[mon-1,:,:,:]
    #convert to kg/s from kg\ms\s
    Jan_loss2 = CO14_loss09.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_loss1 = Jan_loss2.collapsed('time',iris.analysis.MEAN)
    Jan_loss = Jan_loss1.data * 3.15e7 * 30e-3
    #calculating net production and loss 
    Jan_pl = Jan_prod - Jan_loss 
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/Barrier2009/barrier{}.csv'.format(mon-1),delimiter=',', unpack=1)

    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        angle = (180/144)*lati
        for longi in longitudes:
            for h in range(0,85):
                if h >= barrier_height[longi,lati] and angle>=90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi,lati] and angle<90:
                    month_topN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle<90:
                    month_botN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle>=90:
                    month_botS[h,lati,longi]=1
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topN
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topN
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 
    
    Jan_pl_SS = Jan_pl * month_topS
    Jan_pl_SN = Jan_pl * month_topN
    Jan_pl_TS = Jan_pl * month_botS 
    Jan_pl_TN = Jan_pl * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)#/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)#/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)#/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)#/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)#/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)#/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)#/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)#)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
    #finding average over each quadrant for loss
    JanPLSS_avg = np.sum(Jan_pl_SS)#/np.sum(month_topS)
    mean_pl_SS.append(JanPLSS_avg)
    JanPLSN_avg = np.sum(Jan_pl_SN)#/np.sum(month_topN)
    mean_pl_SN.append(JanPLSN_avg)
    JanPLTN_avg = np.sum(Jan_pl_TN)#/np.sum(month_botN)
    mean_pl_TN.append(JanPLTN_avg)
    JanPLTS_avg = np.sum(Jan_pl_TS)#)/np.sum(month_botS)
    mean_pl_TS.append(JanPLTS_avg)
#%% Creating monthly 4 quadrants mask 2010
#finding model level number at which data in tropopause for each month 
print('2010')
for mon in range(1,13):
    print('Month # {} contracts'.format(mon))
    # obtaining production and loss data in Jan 
    Jan_prod2 = CO14_prod10.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_prod1 = Jan_prod2.data * 3.15e7
    Jan_prod = Jan_prod1  * iris.analysis.cartography.area_weights(CO14_prod10)[mon-1,:,:,:]
    #convert to kg/s from kg\ms\s
    Jan_loss2 = CO14_loss10.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_loss1 = Jan_loss2.collapsed('time',iris.analysis.MEAN)
    Jan_loss = Jan_loss1.data * 3.15e7 * 30e-3
    #calculating net production and loss 
    Jan_pl = Jan_prod - Jan_loss 
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/Barrier2010/barrier{}.csv'.format(mon-1),delimiter=',', unpack=1)

    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        angle = (180/144)*lati
        for longi in longitudes:
            for h in range(0,85):
                if h >= barrier_height[longi,lati] and angle>=90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi,lati] and angle<90:
                    month_topN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle<90:
                    month_botN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle>=90:
                    month_botS[h,lati,longi]=1
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topN
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topN
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 
    
    Jan_pl_SS = Jan_pl * month_topS
    Jan_pl_SN = Jan_pl * month_topN
    Jan_pl_TS = Jan_pl * month_botS 
    Jan_pl_TN = Jan_pl * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)#/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)#/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)#/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)#/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)#/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)#/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)#/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)#)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
    #finding average over each quadrant for loss
    JanPLSS_avg = np.sum(Jan_pl_SS)#/np.sum(month_topS)
    mean_pl_SS.append(JanPLSS_avg)
    JanPLSN_avg = np.sum(Jan_pl_SN)#/np.sum(month_topN)
    mean_pl_SN.append(JanPLSN_avg)
    JanPLTN_avg = np.sum(Jan_pl_TN)#/np.sum(month_botN)
    mean_pl_TN.append(JanPLTN_avg)
    JanPLTS_avg = np.sum(Jan_pl_TS)#)/np.sum(month_botS)
    mean_pl_TS.append(JanPLTS_avg)

#%% Creating monthly 4 quadrants mask 2011
#finding model level number at which data in tropopause for each month 
print('2011')
for mon in range(1,13):
    print('Month # {} contracts'.format(mon))
    # obtaining production and loss data in Jan 
    Jan_prod2 = CO14_prod11.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_prod1 = Jan_prod2.data * 3.15e7
    Jan_prod = Jan_prod1  * iris.analysis.cartography.area_weights(CO14_prod10)[mon-1,:,:,:]
    #convert to kg/s from kg\ms\s
    Jan_loss2 = CO14_loss11.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_loss1 = Jan_loss2.collapsed('time',iris.analysis.MEAN)
    Jan_loss = Jan_loss1.data * 3.15e7 * 30e-3
    #calculating net production and loss 
    Jan_pl = Jan_prod - Jan_loss 
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/Barrier2011/barrier{}.csv'.format(mon-1),delimiter=',', unpack=1)

    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        angle = (180/144)*lati
        for longi in longitudes:
            for h in range(0,85):
                if h >= barrier_height[longi,lati] and angle>=90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi,lati] and angle<90:
                    month_topN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle<90:
                    month_botN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle>=90:
                    month_botS[h,lati,longi]=1
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topN
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topN
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 
    
    Jan_pl_SS = Jan_pl * month_topS
    Jan_pl_SN = Jan_pl * month_topN
    Jan_pl_TS = Jan_pl * month_botS 
    Jan_pl_TN = Jan_pl * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)#/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)#/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)#/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)#/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)#/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)#/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)#/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)#)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
    #finding average over each quadrant for loss
    JanPLSS_avg = np.sum(Jan_pl_SS)#/np.sum(month_topS)
    mean_pl_SS.append(JanPLSS_avg)
    JanPLSN_avg = np.sum(Jan_pl_SN)#/np.sum(month_topN)
    mean_pl_SN.append(JanPLSN_avg)
    JanPLTN_avg = np.sum(Jan_pl_TN)#/np.sum(month_botN)
    mean_pl_TN.append(JanPLTN_avg)
    JanPLTS_avg = np.sum(Jan_pl_TS)#)/np.sum(month_botS)
    mean_pl_TS.append(JanPLTS_avg)


#%% Creating monthly 4 quadrants mask 2012
#finding model level number at which data in tropopause for each month 
print('2012')
for mon in range(1,13):
    print('Month # {} contracts'.format(mon))
    # obtaining production and loss data in Jan 
    Jan_prod2 = CO14_prod12.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_prod1 = Jan_prod2.data * 3.15e7
    Jan_prod = Jan_prod1  * iris.analysis.cartography.area_weights(CO14_prod10)[mon-1,:,:,:]
    #convert to kg/s from kg\ms\s
    Jan_loss2 = CO14_loss12.extract(iris.Constraint(time = lambda cell: cell.point.month == mon))
    Jan_loss1 = Jan_loss2.collapsed('time',iris.analysis.MEAN)
    Jan_loss = Jan_loss1.data * 3.15e7 * 30e-3
    #calculating net production and loss 
    Jan_pl = Jan_prod - Jan_loss 
    #converting to kg/s from moles/s
    barrier_height = np.loadtxt('Barrier/Barrier2012/Barrier2012barrier{}.csv'.format(mon-1),delimiter=',', unpack=1)

    #dividing into 4 quadrants 
    month_topS = np.zeros([85,144,192])
    month_topN = np.zeros([85,144,192])
    month_botS = np.zeros([85,144,192])
    month_botN = np.zeros([85,144,192])

    latitudes = np.arange(0,143)
    longitudes = np.arange(0,191)

    for lati in latitudes:
        angle = (180/144)*lati
        for longi in longitudes:
            for h in range(0,85):
                if h >= barrier_height[longi,lati] and angle>=90 :
                    month_topS[h,lati,longi]=1
                if h >= barrier_height[longi,lati] and angle<90:
                    month_topN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle<90:
                    month_botN[h,lati,longi]=1
                if h < barrier_height[longi,lati] and angle>=90:
                    month_botS[h,lati,longi]=1
    #splittin into quadrants 
    Jan_prod_SS = Jan_prod * month_topS
    Jan_prod_SN = Jan_prod * month_topN
    Jan_prod_TS = Jan_prod * month_botS 
    Jan_prod_TN = Jan_prod * month_botN 

    Jan_loss_SS = Jan_loss * month_topS
    Jan_loss_SN = Jan_loss * month_topN
    Jan_loss_TS = Jan_loss * month_botS 
    Jan_loss_TN = Jan_loss * month_botN 
    
    Jan_pl_SS = Jan_pl * month_topS
    Jan_pl_SN = Jan_pl * month_topN
    Jan_pl_TS = Jan_pl * month_botS 
    Jan_pl_TN = Jan_pl * month_botN 

    #finding average over each quadrant for production 
    JanPSS_avg = np.sum(Jan_prod_SS)#/np.sum(month_topS)
    mean_prod_SS.append(JanPSS_avg)
    JanPSN_avg = np.sum(Jan_prod_SN)#/np.sum(month_topN)
    mean_prod_SN.append(JanPSN_avg)
    JanPTN_avg = np.sum(Jan_prod_TN)#/np.sum(month_botN)
    mean_prod_TN.append(JanPTN_avg)
    JanPTS_avg = np.sum(Jan_prod_TS)#/np.sum(month_botS)
    mean_prod_TS.append(JanPTS_avg)

    #finding average over each quadrant for loss
    JanLSS_avg = np.sum(Jan_loss_SS)#/np.sum(month_topS)
    mean_loss_SS.append(JanLSS_avg)
    JanLSN_avg = np.sum(Jan_loss_SN)#/np.sum(month_topN)
    mean_loss_SN.append(JanLSN_avg)
    JanLTN_avg = np.sum(Jan_loss_TN)#/np.sum(month_botN)
    mean_loss_TN.append(JanLTN_avg)
    JanLTS_avg = np.sum(Jan_loss_TS)#)/np.sum(month_botS)
    mean_loss_TS.append(JanLTS_avg)
    
    #finding average over each quadrant for loss
    JanPLSS_avg = np.sum(Jan_pl_SS)#/np.sum(month_topS)
    mean_pl_SS.append(JanPLSS_avg)
    JanPLSN_avg = np.sum(Jan_pl_SN)#/np.sum(month_topN)
    mean_pl_SN.append(JanPLSN_avg)
    JanPLTN_avg = np.sum(Jan_pl_TN)#/np.sum(month_botN)
    mean_pl_TN.append(JanPLTN_avg)
    JanPLTS_avg = np.sum(Jan_pl_TS)#)/np.sum(month_botS)
    mean_pl_TS.append(JanPLTS_avg)

#%% creatiing plots of averages for loss and production
month_numbers = np.arange(1,61)

fig,ax = plt.subplots(2,1, sharex=True)

ax[0].plot(month_numbers,mean_prod_SS, 'x',label='Strat South',color='red')
ax[0].plot(month_numbers,mean_prod_SN, 'x',label='Strat North',color='purple')
ax[0].set_ylim([4.5,7.5])
ax[0].grid()
ax[0].legend(loc ='upper right')

ax[1].plot(month_numbers,mean_prod_TN, 'x',label='Trop North',color='green')
ax[1].plot(month_numbers,mean_prod_TS, 'x',label='Trop South',color='blue')
ax[1].grid()
ax[1].set_ylim([2,3])
ax[1].legend(loc ='upper right')

fig.supxlabel('Month')
fig.supylabel('CO14 Production (kg year${-1}$)')

plt.savefig('Average_Prodall.pdf')
plt.show()

plt.figure()
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].plot(month_numbers,mean_loss_SS,label='Strat South',color='red')
ax[0].plot(month_numbers,mean_loss_SN,label='Strat North',color='purple')
ax[0].grid()
ax[0].set_ylim([1,8])
ax[0].legend(loc ='upper right')

ax[1].plot(month_numbers,mean_loss_TN,label='Trop North',color='green')
ax[1].plot(month_numbers,mean_loss_TS,label='Trop South',color='blue')
ax[1].legend(loc='upper right')
ax[1].set_ylim([1,8])
ax[1].grid()

fig.supxlabel('Month')
fig.supylabel('CO14 Loss (kg year${-1}$)')

plt.savefig('Average_Lossall.pdf')
plt.show()

plt.figure()
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].plot(month_numbers,mean_pl_SN,label='Strat North',color='purple')
ax[0].plot(month_numbers,mean_pl_SS,label='Strat South',color='red')
ax[0].grid()
#ax[0].set_ylim([1,8])
ax[0].legend(loc ='upper right')

ax[1].plot(month_numbers,mean_pl_TN, label='Trop North',color='green')
ax[1].plot(month_numbers,mean_pl_TS, label='Trop South',color='blue')
ax[1].legend(loc='lower right')
#ax[1].set_ylim([1,8])
ax[1].grid()

fig.supxlabel('Month')
fig.supylabel('CO14 Net Prod and Loss (kg year${-1}$)')

plt.savefig('Average_PLall.pdf')
plt.show()

#%%
####
####

#South vs North
####
####
#findining total production for each year in south 
totalprodS08= (sum(mean_prod_SS[0:12])/12)+sum(mean_prod_TS[0:12])/12
totalprodS09= (sum(mean_prod_SS[12:24])/12)+sum(mean_prod_TS[12:24])/12
totalprodS10= (sum(mean_prod_SS[24:36])/12)+sum(mean_prod_TS[24:36])/12
totalprodS11= (sum(mean_prod_SS[36:48])/12)+sum(mean_prod_TS[36:48])/12
totalprodS12= (sum(mean_prod_SS[48:60])/12)+sum(mean_prod_TS[48:60])/12

#findining total production for each year in south 
totalprodN08= (sum(mean_prod_SN[0:12])/12)+sum(mean_prod_TN[0:12])/12
totalprodN09= (sum(mean_prod_SN[12:24])/12)+sum(mean_prod_TN[12:24])/12
totalprodN10= (sum(mean_prod_SN[24:36])/12)+sum(mean_prod_TN[24:36])/12
totalprodN11= (sum(mean_prod_SN[36:48])/12)+sum(mean_prod_TN[36:48])/12
totalprodN12= (sum(mean_prod_SN[48:60])/12)+sum(mean_prod_TN[48:60])/12

#findining total loss for each year in south
totallossS08= (sum(mean_loss_SS[0:12])/12)+sum(mean_loss_TS[0:12])/12
totallossS09= (sum(mean_loss_SS[12:24])/12)+sum(mean_loss_TS[12:24])/12
totallossS10= (sum(mean_loss_SS[24:36])/12)+sum(mean_loss_TS[24:36])/12
totallossS11= (sum(mean_loss_SS[36:48])/12)+sum(mean_loss_TS[36:48])/12
totallossS12= (sum(mean_loss_SS[48:60])/12)+sum(mean_loss_TS[48:60])/12

#findining total loss for each year in nort2
totallossN08= (sum(mean_loss_SN[0:12])/12)+sum(mean_loss_TN[0:12])/12
totallossN09= (sum(mean_loss_SN[12:24])/12)+sum(mean_loss_TN[12:24])/12
totallossN10= (sum(mean_loss_SN[24:36])/12)+sum(mean_loss_TN[24:36])/12
totallossN11= (sum(mean_loss_SN[36:48])/12)+sum(mean_loss_TN[36:48])/12
totallossN12= (sum(mean_loss_SN[48:60])/12)+sum(mean_loss_TN[48:60])/12

#####
#####
#Trop vs Strat
#####
##### 
#findining total production for each year in TROP
totalprodT08= (sum(mean_prod_TN[0:12])/12)+sum(mean_prod_TS[0:12])/12
totalprodT09= (sum(mean_prod_TN[12:24])/12)+sum(mean_prod_TS[12:24])/12
totalprodT10= (sum(mean_prod_TN[24:36])/12)+sum(mean_prod_TS[24:36])/12
totalprodT11= (sum(mean_prod_TN[36:48])/12)+sum(mean_prod_TS[36:48])/12
totalprodT12= (sum(mean_prod_TN[48:60])/12)+sum(mean_prod_TS[48:60])/12

#findining total production for each year in STRAT
totalprodST08= (sum(mean_prod_SN[0:12]))/12+sum(mean_prod_SS[0:12])/12
totalprodST09= (sum(mean_prod_SN[12:24])/12)+sum(mean_prod_SS[12:24])/12
totalprodST10= (sum(mean_prod_SN[24:36])/12)+sum(mean_prod_SS[24:36])/12
totalprodST11= (sum(mean_prod_SN[36:48])/12)+sum(mean_prod_SS[36:48])/12
totalprodST12= (sum(mean_prod_SN[48:60])/12)+sum(mean_prod_SS[48:60])/12

#findining total loss for each year in TROP
totallossT08= (sum(mean_loss_TN[0:12])/12) + sum(mean_loss_TS[0:12])/12
totallossT09= (sum(mean_loss_TN[12:24])/12) +sum(mean_loss_TS[12:24])/12
totallossT10= (sum(mean_loss_TN[24:36])/12)+sum(mean_loss_TS[24:36])/12
totallossT11= (sum(mean_loss_TN[36:48])/12)+sum(mean_loss_TS[36:48])/12
totallossT12= (sum(mean_loss_TN[48:60])/12)+sum(mean_loss_TS[48:60])/12

#findining total loss for each year in STRAT
totallossST08= (sum(mean_loss_SN[0:12])/12)+sum(mean_loss_SS[0:12])/12
totallossST09= (sum(mean_loss_SN[12:24])/12)+sum(mean_loss_SS[12:24])/12
totallossST10= (sum(mean_loss_SN[24:36])/12)+sum(mean_loss_SS[24:36])/12
totallossST11= (sum(mean_loss_SN[36:48])/12)+sum(mean_loss_SS[36:48])/12
totallossST12= (sum(mean_loss_SN[48:60])/12)+sum(mean_loss_SS[48:60])/12


#FINDING TOTAL PRODUCTION AND LOSS FOR EACH YEAR 
total_prod_annual =[]
total_loss_annual =[]

prod_08 = totalprodS08 + totalprodN08 
prod_09 = totalprodS09 + totalprodN09
prod_10 = totalprodS10 + totalprodN10 
prod_11 = totalprodS11 + totalprodN11 
prod_12 = totalprodS12 + totalprodN12 

loss_08 = totallossS08 + totallossN08
loss_09 = totallossS09 + totallossN09
loss_10 = totallossS10 + totallossN10 
loss_11 = totallossS11 + totallossN11 
loss_12 = totallossS12 + totallossN12 

total_prod_annual =[prod_08,prod_09,prod_10,prod_11,prod_12]
total_loss_annual =[loss_08,loss_09,loss_10,loss_11,loss_12]


plt.plot([2008,2009,2010,2011,2012],[netprod08,netprod09,netprod10,netprod11,netprod12],'x',label='Exp Prod',color='red')
plt.plot([2008,2009,2010,2011,2012],[netloss08,netloss09,netloss10,netloss11,netloss12],'x',label='Exp Loss',color='blue')
plt.plot([2008,2009,2010,2011,2012],np.array(total_prod_annual),'x',label='Mes Prod',color='orange')
plt.plot([2008,2009,2010,2011,2012],np.array(total_loss_annual),'x',label='Mes Loss',color='purple')

plt.xticks([2008,2009,2010,2011,2012])
plt.xlabel('Year')
plt.ylabel('Loss or Production (kg year ${-1}$)')
plt.title('Annual Loss and Production')
plt.grid()
plt.legend()
plt.show()

prod_trop =  [totalprodT08,totalprodT09,totalprodT10,totalprodT11,totalprodT12]
loss_trop =  [totallossT08,totallossT09,totallossT10,totallossT11,totallossT12]

prod_strat = [totalprodST08,totalprodST09,totalprodST10,totalprodST11,totalprodST12]
loss_strat = [totallossST08,totallossST09,totallossST10,totallossST11,totallossST12]


prod_south = [totalprodS08,totalprodS09,totalprodS10,totalprodS11,totalprodS12]
loss_south = [totallossS08,totallossS09,totallossS10,totallossS11,totallossS12]

prod_north = [totalprodN08,totalprodN09,totalprodN10,totalprodN11,totalprodN12]
loss_north = [totallossN08,totallossN09,totallossN10,totallossN11,totallossN12]

plt.figure()
fig,ax = plt.subplots(2,1, sharex=True)
ax[0].plot([2008,2009,2010,2011,2012],prod_strat,'x',label='S',color='blue')
ax[0].plot([2008,2009,2010,2011,2012],prod_trop,'x',label='T',color='red')
ax[0].set_ylabel('Production (kg year ${-1}$)')
ax[0].grid()
ax[0].set_ylim([4,14])
ax[0].legend()

ax[1].plot([2008,2009,2010,2011,2012],loss_strat,'x',label='S',color='blue')
ax[1].plot([2008,2009,2010,2011,2012],loss_trop,'x',label='T',color='red')
ax[1].set_ylabel('Loss (kg year ${-1}$)')
ax[1].grid()
ax[1].set_ylim([7,10])
ax[1].legend()
fig.supxlabel('Year')
plt.xticks([2008,2009,2010,2011,2012])
plt.show()

plt.figure()
fig,ax = plt.subplots(2,1, sharex=True)

ax[0].plot([2008,2009,2010,2011,2012],prod_south,'x',label='South',color='blue')
ax[0].plot([2008,2009,2010,2011,2012],prod_north,'x',label='North',color='red')
ax[0].set_ylabel('Production (kg year ${-1}$)')
ax[0].set_ylim([7,10])
ax[0].legend()
ax[0].grid()

ax[1].plot([2008,2009,2010,2011,2012],loss_south,'x',label='South',color='blue')
ax[1].plot([2008,2009,2010,2011,2012],loss_north,'x',label='North',color='red')
ax[1].set_ylabel('Loss (kg year ${-1}$)')
ax[1].set_ylim([7,10])
ax[1].legend()
ax[1].grid()


plt.xticks([2008,2009,2010,2011,2012])
fig.supxlabel('Year')
plt.show()

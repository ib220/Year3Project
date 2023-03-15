#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 13:31:18 2023

@author: maryamfatima
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd


allweights = pd.read_csv('/Users/maryamfatima/Desktop/weights.csv')


PSN=allweights.iloc[0]

PSS=allweights.iloc[1]
PSTN=allweights.iloc[2]
PSTS=allweights.iloc[3]

PTN=allweights.iloc[4]
PTS=allweights.iloc[5]
PTTN=allweights.iloc[6]
PTTS=allweights.iloc[7]



LSN=allweights.iloc[8]
LSS=allweights.iloc[9]
LSTN=allweights.iloc[10]
LSTS=allweights.iloc[11]

LTN=allweights.iloc[12]
LTS=allweights.iloc[13]
LTTN=allweights.iloc[14]
LTTS=allweights.iloc[15]



PLSN=allweights.iloc[16]
PLSS=allweights.iloc[17]
PLSTN=allweights.iloc[18]
PLSTS=allweights.iloc[19]

PLTN=allweights.iloc[20]
PLTS=allweights.iloc[21]
PLTTN=allweights.iloc[22]
PLTTS=allweights.iloc[23]
    
    
#%%


month_numbers = np.arange(1,61)


fig,ax = plt.subplots(2,1, sharex=True)

plt.title('Loss')

ax[0].plot(month_numbers,LSTS,label="stratosphere tropic south",color='red')
ax[0].plot(month_numbers,LSTN,label="stratosphere tropic north",color='purple')
ax[0].plot(month_numbers,LSS,label="stratosphere extratropic south",color='orange')
ax[0].plot(month_numbers,LSN,label="stratosphere extratropic north",color='violet')

ax[0].grid()
ax[0].set_ylim([0,7])
ax[0].legend(fontsize=10)

ax[1].plot(month_numbers,LTTS,label="trposphere tropic south",color='black')
ax[1].plot(month_numbers,LTTN,label="troposphere tropic north",color='teal')
ax[1].plot(month_numbers,LTS,label="troposphere extratropic south",color='green')
ax[1].plot(month_numbers,LTN,label="troposphere extratropic north",color='blue')


ax[1].grid()
ax[1].set_ylim([0,7])
ax[1].legend(fontsize=10)

plt.savefig('/Users/maryamfatima/Desktop/loss_sectors.png')



#%%
PL_total_strat_north=PLSN+PLSTN
PL_total_strat_south=PLSS+PLSTS
PL_total_trop_north=PLTN+PLTTN
PL_total_trop_south=PLTS+PLTTS


L_total_strat_north=LSN+LSTN
L_total_strat_south=LSS+LSTS
L_total_trop_north=LTN+LTTN
L_total_trop_south=LTS+LTTS

P_total_strat_north=PSN+PSTN
P_total_strat_south=PSS+PSTS
P_total_trop_north=PTN+PTTN
P_total_trop_south=PTS+PTTS

PL_ratio_north=PL_total_strat_north/PL_total_trop_north
PL_ratio_south=PL_total_strat_south/PL_total_trop_south


L_ratio_north=L_total_strat_north/L_total_trop_north
L_ratio_south=L_total_strat_south/L_total_trop_south

P_ratio_north=P_total_strat_north/P_total_trop_north
P_ratio_south=P_total_strat_south/P_total_trop_south




PL_ratio_TN=PLSTN/PLTTN
PL_ratio_TS=PLSTS/PLTTS
PL_ratio_N=PLSN/PLTN
PL_ratio_S=PLSS/PLTS

L_ratio_TN=LSTN/LTTN
L_ratio_TS=LSTS/LTTS
L_ratio_N=LSN/LTN
L_ratio_S=LSS/LTS

P_ratio_TN=PSTN/PTTN
P_ratio_TS=PSTS/PTTS
P_ratio_N=PSN/PTN
P_ratio_S=PSS/PTS


#%%

fig,ax = plt.subplots(2,1, sharex=True)

plt.title('Loss Ratios between Stratosphere and Troposphere')

ax[0].plot(month_numbers,L_ratio_N,label="Extratropic North",color='red')
ax[0].plot(month_numbers,L_ratio_TN,label="Tropic North",color='brown')

ax[0].grid()
#ax[0].set_ylim([0,2])
ax[0].legend(fontsize=10)

ax[1].plot(month_numbers,L_ratio_S,label="Extratropic South",color='blue')
ax[1].plot(month_numbers,L_ratio_TS,label="Tropic South",color='green')


ax[1].grid()
#ax[1].set_ylim([0,2])
ax[1].legend(fontsize=10)

plt.savefig('/Users/maryamfatima/Desktop/loss_ratios.png')

#%%

av_L_rat_N=np.mean(L_ratio_N)
av_L_rat_S=np.mean(L_ratio_S)
av_L_rat_TN=np.mean(L_ratio_TN)
av_L_rat_TS=np.mean(L_ratio_TS)

av_P_rat_S=np.mean(P_ratio_S)
av_P_rat_N=np.mean(P_ratio_N)
av_P_rat_TS=np.mean(P_ratio_TS)
av_P_rat_TN=np.mean(P_ratio_TN)

av_PL_rat_N=np.mean(PL_ratio_N)
av_PL_rat_S=np.mean(PL_ratio_S)
av_PL_rat_TN=np.mean(PL_ratio_TN)
av_PL_rat_TS=np.mean(PL_ratio_TS)


print(av_L_rat_N)
print(av_L_rat_S)
print(av_L_rat_TN)
print(av_L_rat_TS)

print("\n")

print(av_P_rat_S)
print(av_P_rat_N)
print(av_P_rat_TS)
print(av_P_rat_TN)

print("\n")

print(av_PL_rat_N)
print(av_PL_rat_S)
print(av_PL_rat_TN)
print(av_PL_rat_TS)





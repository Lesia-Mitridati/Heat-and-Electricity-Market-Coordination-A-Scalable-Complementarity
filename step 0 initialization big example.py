# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 08:57:30 2016

@author: lemitri
"""

#%%
            
import os
import pandas as pd
import scipy.stats as sp
#import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
sb.set_style('ticks')

import gurobipy as gb
import itertools as it

import numpy as np

os.chdir("C:/Users/lemitri/Documents/PhD/Sequential and Hierarchical Heat and Electricity Markets/python/EJOR article/big example")

from sklearn.cluster import KMeans

W_range=[50,75,100,125,150,175,200,225,250,275,300]
 
# indexes
    
T=24
G=7 #elec generator
W=6 #wind
CH=4 #chp
HO=1 #heat only
HS=3  #3 #storage
ES=0 # elec storage
HP=3 #heat pump

S_load=15
S_wind=15
S_supply=15
S_all=S_load*S_wind*S_supply


day = ['d{0:02d}'.format(t) for t in range(31)]
time = ['t{0:02d}'.format(t) for t in range(T)]
time_list=np.arange(T)
gen=['G{0}'.format(g+1) for g in range(G)]
heat_storage=['HS{0}'.format(s+1) for s in range(HS)]
elec_storage=['ES{0}'.format(s+1) for s in range(ES)]
heat_pump = ['HP{0}'.format(h+1) for h in range(HP)]
wind=['W{0}'.format(w+1) for w in range(W)]
heat_only = ['HO{0}'.format(ho+1) for ho in range(HO)]
#heat_exchanger_station = ['HES{0}'.format(h+1) for h in range(HES)]
CHP_sorted = {'ex':['CHP1','CHP2','CHP3','CHP4'],'bp':[]} # CHPs indexes sorted by type: extraction or backpressure
CHP = list(it.chain.from_iterable(CHP_sorted.values()))

scenario_load=['S_load{0}'.format(s+1) for s in range(S_load)]
scenario_wind=['S_wind{0}'.format(s+1) for s in range(S_wind)]
scenario_supply=['S_supply{0}'.format(s+1) for s in range(S_supply)]
scenario=['S{0}'.format(s+1) for s in range(S_all)]
scenario_list=[[scenario_load[s1],scenario_wind[s2],scenario_supply[s3]] for s1 in range(S_load) for s2 in range(S_wind) for s3 in range(S_supply)]
scenario_dic={s:x for (s,x) in zip(scenario,scenario_list)}    
          
#%% loads
           
heat=pd.read_csv("total heat load dec 2017.csv",sep=";",decimal=".")

heat_load_dec={}
for d in range(31):
    for t in range(T):
        heat_load_dec[day[d],time[t]]=heat['TOTAL'][24*d+t]
        
heat_load_average = {t:(heat_load_dec[day[0],t]+heat_load_dec[day[2],t]+heat_load_dec[day[3],t])/3 for t in time}


heat_load_norm = {t:(heat_load_average[t]-min([heat_load_average[t] for t in time]))/(max([heat_load_average[t] for t in time])-min([heat_load_average[t] for t in time])) for t in time}

heat_load = {}

min1= 450 # 1500/2
max1= 1000 # 2500/2
for t in time:
    heat_load[t]=heat_load_norm[t]*(max1-min1)+min1
      
#%%  
      
elec_load_IEEE = {'t00':1775.835,'t01':1669.815,'t02':1590.3,'t03':1563.795,'t04':1563.795,'t05':1590.3,'t06':1961.37,'t07':2279.43 ,'t08':2517.975,'t09':2544.48,'t10':2544.48,'t11':2517.975,'t12':2517.975,'t13':2517.975,'t14':2464.965,'t15':2464.965,'t16':2623.995,'t17':2650.5,'t18':2650.5,'t19':2544.48,'t20':2411.955,'t21':2199.915,'t22':1934.865,'t23':1669.815}
elec_load_IEEE_min=min([elec_load_IEEE[t] for t in time])
elec_load_IEEE_max=max([elec_load_IEEE[t] for t in time])

elec=pd.read_csv("total elec load jan 2018.csv",sep=";",decimal=",")
#heat_load = {'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}

elec_load_jan={}
for d in range(31):
    for t in range(T):
        elec_load_jan[day[d],time[t]]=elec['Total'][24*d+t]
 
elec_load_average = {t:sum(elec_load_jan[d,t] for d in day)/31 for t in time}

elec_load_scenario_0 = {}
for s in scenario_load:
    dd=np.random.randint(0,31)    
    for t in time:
        elec_load_scenario_0[s,t]= elec_load_jan[day[dd],t]

elec_load_scenario_norm = {(s,t):(elec_load_scenario_0[s,t]-min([elec_load_scenario_0[s,t] for t in time]))/(max([elec_load_scenario_0[s,t] for t in time])-min([elec_load_scenario_0[s,t] for t in time])) for s in scenario_load for t in time}


#%%

elec_load_scenario = {}
min1 = {}
max1 = {}

for s in scenario_load:
#    min1[s]=sp.norm(elec_load_IEEE_min,elec_load_IEEE_min*0.05).rvs()
#    max1[s]=sp.norm(elec_load_IEEE_max,elec_load_IEEE_min*0.05).rvs()
    min1[s]= 600 # 1500/2
    max1[s]= 1300 # 2500/2
    for t in time:
        elec_load_scenario[s,t]=elec_load_scenario_norm[s,t]*(max1[s]-min1[s])+min1[s]


#%% cost

alpha_m = {}
alpha_v = {}

alpha_m['HO1']=100 # TRY WITH 60???


alpha_m['G5'] = 6.02 
alpha_m['G6'] = 5.47

alpha_m['CHP1']=5 # elec = 30
alpha_m['CHP4']=7.5
alpha_m['CHP2']=10 # elec = 11.4
alpha_m['CHP3']=12.5 # elec = 19

alpha_m['G1'] = 13.32
alpha_m['G2'] = 13.32 


alpha_m['G3'] = 20.7 
alpha_m['G4'] = 26.11 


alpha_m['G7'] = 50.0


alpha_m['W1']=0.0000001
alpha_m['W2']=0.0000001
alpha_m['W3']=0.0000001
alpha_m['W4']=0.0000001
alpha_m['W5']=0.0000001
alpha_m['W6']=0.0000001

for g in CHP+heat_only+heat_pump:
    alpha_v[g] = 0.0

for g in gen:
    alpha_v[g] = alpha_m[g]*0.1

alpha_v['G7'] = 0.0
alpha_v['W1']=0.0
alpha_v['W2']=0.0
alpha_v['W3']=0.0
alpha_v['W4']=0.0
alpha_v['W5']=0.0
alpha_v['W6']=0.0

alpha_scenario = {}

for s in scenario_supply:
    for t in time:
        for g in CHP+heat_only+gen+wind:
            alpha_scenario[s,g,t]=sp.norm(alpha_m[g],alpha_v[g]).rvs()
                  
#%% technical characteristics
          
heat_maxprod = {'CHP1': 300,'CHP2': 300,'CHP3': 300,'CHP4': 200,'HP1':250,'HP2':250,'HP3':250,'HO1':1000} #HP 150

rho_elec = {'CHP1': 2.1,'CHP2': 2.1,'CHP3': 2.1,'CHP4':2.4} # efficiency of the CHP for electricity production
rho_heat = {'CHP1': 0.25,'CHP2': 0.21,'CHP3': 0.25,'CHP4':0.21,'HO1':1} # efficiency of the CHP for heat production

r_min = {'CHP1': 0.5,'CHP2': 0.5,'CHP3': 0.5,'CHP4':0.5} # elec/heat ratio (flexible in the case of extraction units) 
CHP_maxprod = {'CHP1':600,'CHP2':600,'CHP3':600,'CHP4':600}

COP_heat_pump={'HP1':2.8,'HP2':3.1,'HP3':2.5}
COP = COP_heat_pump

storage_loss={h:1.5 for h in heat_storage+elec_storage}
storage_init={h:100 for h in heat_storage+elec_storage}
storage_rho_plus={h:1.1 for h in heat_storage+elec_storage} # >=1
storage_rho_moins={h:0.9 for h in heat_storage+elec_storage} # <=1
storage_maxcapacity={h:150 for h in heat_storage+elec_storage}
storage_maxprod = {h:50 for h in heat_storage+elec_storage}

alpha_heat_m = {h:rho_heat[h]*alpha_m[h] for h in CHP}
alpha_elec_m = {h:rho_elec[h]*alpha_m[h] for h in CHP}

elec_maxprod = {} # known
for g in CHP_sorted['ex']:
    elec_maxprod[g] = CHP_maxprod[g]/rho_elec[g]
for g in CHP_sorted['bp']:
    elec_maxprod[g] = heat_maxprod[g]*r_min[g]   
    
for g in wind: 
    elec_maxprod[g] = 1000 # same ratio = 500/2
    
elec_maxprod['G1'] = 152/2
elec_maxprod['G2'] = 152/2
elec_maxprod['G3'] = 350/2
elec_maxprod['G4'] = 60/2
elec_maxprod['G5'] = 400/2
elec_maxprod['G6'] = 400/2
elec_maxprod['G7'] = 1000


#%% wind scenarios

wind_file = pd.Panel({'W' + str(i+1):pd.read_csv('wind {0}.out'.format(i+1), index_col=0, skip_footer=True) for i in range(W)})
           
wind_scenario = {}
for w in wind:
    for t in range(T):
        for s in range(S_wind):
            wind_scenario[scenario_wind[s],w,time[t]] = wind_file[w,t+1,'V{0}'.format(s+1)]

#%% scenario reduction 

S_kmeans_0 = 1
S_all_kmeans_0 = 1

scenario_load_kmeans_0=['S_load{0}'.format(s+1) for s in range(S_kmeans_0)]
scenario_wind_kmeans_0=['S_wind{0}'.format(s+1) for s in range(S_kmeans_0)]
scenario_supply_kmeans_0=['S_supply{0}'.format(s+1) for s in range(S_kmeans_0)]

scenario_kmeans_0=['S{0}'.format(s+1) for s in range(S_all_kmeans_0)]
scenario_list_kmeans_0=[[scenario_load[s1],scenario_wind[s2],scenario_supply[s3]] for s1 in range(S_kmeans_0) for s2 in range(S_kmeans_0) for s3 in range(S_kmeans_0)]
scenario_dic_kmeans_0={s:x for (s,x) in zip(scenario_kmeans_0,scenario_list_kmeans_0)}    


data_wind = [[wind_scenario[s,'W1',t] for t in time]+[wind_scenario[s,'W2',t] for t in time]+[wind_scenario[s,'W3',t] for t in time]+[wind_scenario[s,'W4',t] for t in time]+[wind_scenario[s,'W5',t] for t in time]+[wind_scenario[s,'W6',t] for t in time] for s in scenario_wind]
data_supply = [[alpha_scenario[s,'G1',t] for t in time]+[alpha_scenario[s,'G2',t] for t in time]+[alpha_scenario[s,'G3',t] for t in time]+[alpha_scenario[s,'G4',t] for t in time]+[alpha_scenario[s,'G5',t] for t in time]+[alpha_scenario[s,'G6',t] for t in time]+[alpha_scenario[s,'G7',t] for t in time] for s in scenario_supply]
data_load = [[elec_load_scenario[s,t] for t in time] for s in scenario_load]

elec_load_scenario_kmeans_0 = {}
alpha_scenario_kmeans_0 = {}
wind_scenario_kmeans_0 = {}

kmeans_0_load = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_load)
kmeans_0_wind = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_wind)
kmeans_0_supply = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_supply)

for s in range(S_kmeans_0): 
    for t in range(T):
        elec_load_scenario_kmeans_0[scenario_load_kmeans_0[s],time[t]]=kmeans_0_load.cluster_centers_[s][t]


for s in range(S_kmeans_0): 
    for t in range(T):
        for g in range(G):
            alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],gen[g],time[t]]=kmeans_0_supply.cluster_centers_[s][t+24*g] 
            
        for h in CHP+heat_only+wind:
            alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario[scenario_supply_kmeans_0[s],h,time[t]] 

for s in range(S_kmeans_0):
    for w in range(W): 
        for t in range(T):                
            wind_scenario_kmeans_0[scenario_wind_kmeans_0[s],wind[w],time[t]]=kmeans_0_wind.cluster_centers_[s][t+24*w]
            
            
#%%

S_kmeans = 6
S_all_kmeans = S_kmeans**3

scenario_load_kmeans=['S_load{0}'.format(s+1) for s in range(S_kmeans)]
scenario_wind_kmeans=['S_wind{0}'.format(s+1) for s in range(S_kmeans)]
scenario_supply_kmeans=['S_supply{0}'.format(s+1) for s in range(S_kmeans)]

scenario_kmeans=['S{0}'.format(s+1) for s in range(S_all_kmeans)]
scenario_list_kmeans=[[scenario_load[s1],scenario_wind[s2],scenario_supply[s3]] for s1 in range(S_kmeans) for s2 in range(S_kmeans) for s3 in range(S_kmeans)]
scenario_dic_kmeans={s:x for (s,x) in zip(scenario_kmeans,scenario_list_kmeans)}    


data_wind = [[wind_scenario[s,'W1',t] for t in time]+[wind_scenario[s,'W2',t] for t in time]+[wind_scenario[s,'W3',t] for t in time]+[wind_scenario[s,'W4',t] for t in time]+[wind_scenario[s,'W5',t] for t in time]+[wind_scenario[s,'W6',t] for t in time] for s in scenario_wind]
data_supply = [[alpha_scenario[s,'G1',t] for t in time]+[alpha_scenario[s,'G2',t] for t in time]+[alpha_scenario[s,'G3',t] for t in time]+[alpha_scenario[s,'G4',t] for t in time]+[alpha_scenario[s,'G5',t] for t in time]+[alpha_scenario[s,'G6',t] for t in time]+[alpha_scenario[s,'G7',t] for t in time] for s in scenario_supply]
data_load = [[elec_load_scenario[s,t] for t in time] for s in scenario_load]

elec_load_scenario_kmeans = {}
alpha_scenario_kmeans = {}
wind_scenario_kmeans = {}

kmeans_load = KMeans(n_clusters=S_kmeans, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_load)
kmeans_wind = KMeans(n_clusters=S_kmeans, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_wind)
kmeans_supply = KMeans(n_clusters=S_kmeans, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_supply)

for s in range(S_kmeans): 
    for t in range(T):
        elec_load_scenario_kmeans[scenario_load_kmeans[s],time[t]]=kmeans_load.cluster_centers_[s][t]


for s in range(S_kmeans): 
    for t in range(T):
        for g in range(G):
            alpha_scenario_kmeans[scenario_supply_kmeans[s],gen[g],time[t]]=kmeans_supply.cluster_centers_[s][t+24*g] 
            
        for h in CHP+heat_only+wind:
            alpha_scenario_kmeans[scenario_supply_kmeans[s],h,time[t]]=alpha_scenario[scenario_supply_kmeans[s],h,time[t]] 

for s in range(S_kmeans):
    for w in range(W): 
        for t in range(T):                
            wind_scenario_kmeans[scenario_wind_kmeans[s],wind[w],time[t]]=kmeans_wind.cluster_centers_[s][t+24*w]

alpha_up_kmeans = {}
alpha_down_kmeans = {}

for s in range(S_kmeans): 
    for t in range(T):
        for h in CHP:
            alpha_up_kmeans[scenario_supply_kmeans[s],h,time[t]]=alpha_scenario_kmeans[scenario_supply_kmeans[s],h,time[t]]*rho_heat[h]*1.1
            alpha_down_kmeans[scenario_supply_kmeans[s],h,time[t]]=alpha_scenario_kmeans[scenario_supply_kmeans[s],h,time[t]]*rho_heat[h]*0.9

        for h in heat_only:
            alpha_up_kmeans[scenario_supply_kmeans[s],h,time[t]]=alpha_scenario_kmeans[scenario_supply_kmeans[s],h,time[t]]*1.1
            alpha_down_kmeans[scenario_supply_kmeans[s],h,time[t]]=alpha_scenario_kmeans[scenario_supply_kmeans[s],h,time[t]]*0.9
   
#%%

S_kmeans_0 = 1
S_all_kmeans_0 = S_kmeans_0**3

scenario_load_kmeans_0=['S_load{0}'.format(s+1) for s in range(S_kmeans_0)]
scenario_wind_kmeans_0=['S_wind{0}'.format(s+1) for s in range(S_kmeans_0)]
scenario_supply_kmeans_0=['S_supply{0}'.format(s+1) for s in range(S_kmeans_0)]

scenario_kmeans_0=['S{0}'.format(s+1) for s in range(S_all_kmeans_0)]
scenario_list_kmeans_0=[[scenario_load[s1],scenario_wind[s2],scenario_supply[s3]] for s1 in range(S_kmeans_0) for s2 in range(S_kmeans_0) for s3 in range(S_kmeans_0)]
scenario_dic_kmeans_0={s:x for (s,x) in zip(scenario_kmeans_0,scenario_list_kmeans_0)}    


data_wind = [[wind_scenario[s,'W1',t] for t in time]+[wind_scenario[s,'W2',t] for t in time]+[wind_scenario[s,'W3',t] for t in time]+[wind_scenario[s,'W4',t] for t in time]+[wind_scenario[s,'W5',t] for t in time]+[wind_scenario[s,'W6',t] for t in time] for s in scenario_wind]
data_supply = [[alpha_scenario[s,'G1',t] for t in time]+[alpha_scenario[s,'G2',t] for t in time]+[alpha_scenario[s,'G3',t] for t in time]+[alpha_scenario[s,'G4',t] for t in time]+[alpha_scenario[s,'G5',t] for t in time]+[alpha_scenario[s,'G6',t] for t in time]+[alpha_scenario[s,'G7',t] for t in time] for s in scenario_supply]
data_load = [[elec_load_scenario[s,t] for t in time] for s in scenario_load]

elec_load_scenario_kmeans_0 = {}
alpha_scenario_kmeans_0 = {}
wind_scenario_kmeans_0 = {}

kmeans_load = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_load)
kmeans_wind = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_wind)
kmeans_supply = KMeans(n_clusters=S_kmeans_0, init='k-means++', n_init=10, max_iter=300, tol=0.0001,verbose=0, random_state=None, copy_x=True).fit(data_supply)

for s in range(S_kmeans_0): 
    for t in range(T):
        elec_load_scenario_kmeans_0[scenario_load_kmeans_0[s],time[t]]=kmeans_load.cluster_centers_[s][t]


for s in range(S_kmeans_0): 
    for t in range(T):
        for g in range(G):
            alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],gen[g],time[t]]=kmeans_supply.cluster_centers_[s][t+24*g] 
            
        for h in CHP+heat_only+wind:
            alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario[scenario_supply_kmeans_0[s],h,time[t]] 

for s in range(S_kmeans_0):
    for w in range(W): 
        for t in range(T):                
            wind_scenario_kmeans_0[scenario_wind_kmeans_0[s],wind[w],time[t]]=kmeans_wind.cluster_centers_[s][t+24*w]

alpha_up_kmeans_0 = {}
alpha_down_kmeans_0 = {}

for s in range(S_kmeans_0): 
    for t in range(T):
        for h in CHP:
            alpha_up_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]*rho_heat[h]*1.1
            alpha_down_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]*rho_heat[h]*0.9

        for h in heat_only:
            alpha_up_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]*1.1
            alpha_down_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]=alpha_scenario_kmeans_0[scenario_supply_kmeans_0[s],h,time[t]]*0.9

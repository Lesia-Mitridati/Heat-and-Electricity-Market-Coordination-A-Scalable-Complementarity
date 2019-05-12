# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:26:59 2016

@author: lemitri
"""

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

os.chdir("C:/Users/lemitri/Documents/phd/Sequential and Hierarchical Heat and Electricity Markets/DATA")


#%% heat marginal costs calculated as a function of spot prices estimates! 

alpha_heat_sequential = {}  

for W_max in W_range:
      
    for h in CHP_sorted['bp']:
        for t in time:
            alpha_heat_sequential[W_max,h,t] = alpha_scenario_kmeans['S_supply1',h,t]*(rho_heat[h] + rho_elec[h]*r_min[h])-spot_price_hierarchical_average_0[W_max,t]*r_min[h]
    
    
    for h in CHP_sorted['ex']:
        for t in time:
            p_estimate=spot_price_hierarchical_average_0[W_max,t]
            if p_estimate < alpha_scenario_kmeans['S_supply1',h,t]*rho_elec[h]:
                alpha_heat_sequential[W_max,h,t] = alpha_scenario_kmeans['S_supply1',h,t]*(rho_heat[h] + rho_elec[h]*r_min[h])-p_estimate*r_min[h]
    
            else:
                alpha_heat_sequential[W_max,h,t] =  p_estimate*rho_heat[h]/rho_elec[h]
    
    for h in heat_pump:  
        for t in time:
            alpha_heat_sequential[W_max,h,t] =  spot_price_hierarchical_average_0[W_max,t]/COP[h]
            
    for h in heat_only:  
        for t in time:
            alpha_heat_sequential[W_max,h,t] = alpha_scenario_kmeans['S_supply1',h,t]


#%% heat dispatch

class expando(object):
    '''
        A small class which can have attributes set
    '''
    pass

class sequential_heat:
    def __init__(self,W_max):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data(W_max)
        self._build_model()
    
    def optimize(self):
        self.model.optimize()
    
    def _load_data(self,W_max):
        
        #indexes
        self.data.time = time
        self.data.time_list=time_list
        self.data.heat_storage=heat_storage
        self.data.heat_pump =heat_pump
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP

        
        # LOADS
        self.data.heat_load = heat_load
        self.data.heat_maxprod = {(g,t):heat_maxprod[g] for g in CHP+heat_only+heat_pump for t in time} #BIDS
        self.data.storage_discharge_eff= storage_rho_plus
        self.data.storage_charge_eff= storage_rho_moins
        self.data.storage_maxcapacity= storage_maxcapacity
        self.data.storage_maxprod= storage_maxprod
        self.data.storage_loss= storage_loss
        self.data.storage_energy_init= storage_init
        
        # Cost parameters
        self.data.alpha = {(g,t):alpha_heat_sequential[W_max,g,t] for g in CHP+heat_only+heat_pump for t in time} #BIDS
   
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
        self.model.update()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.data.time
#        node=self.data.node
#        line=self.data.line
#        pipe=self.data.pipe
        #gen=self.data.gen
        heat_storage=self.data.heat_storage
        #elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        #wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP
        m = self.model
        
        #heat market optimization variables

        self.variables.Q = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only+heat_pump:
                self.variables.Q[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h,t],name='Q({0},{1})'.format(h,t))

        self.variables.storage_discharge = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_discharge[h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage discharge({0},{1})'.format(h,t))
                    
        self.variables.storage_charge = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_charge[h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage charge({0},{1})'.format(h,t))
                    

        self.variables.storage_energy = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                self.variables.storage_energy[h,t] = m.addVar(lb=0,ub=self.data.storage_maxcapacity[h],name='storage energy({0},{1})'.format(h,t))                
 
        m.update()
        
    def _build_objective(self): # building the objective function for the heat maret clearing
        
        #indexes shortcuts 
        time = self.data.time
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted    
        m = self.model
             
        m.setObjective(
            gb.quicksum(self.variables.Q[h,t]*self.data.alpha[h,t] for h in CHP+heat_only+heat_pump for t in time),
            gb.GRB.MINIMIZE)
            
        
    def _build_constraints(self):
        

        #indexes shortcuts 
        time = self.data.time
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP
        m = self.model
        
        ##1) heat balance equation
        
        self.constraints.heat_balance = {}
        
        for t in time:
                
                self.constraints.heat_balance[t] = m.addConstr(
                    gb.quicksum(self.variables.storage_discharge[i,t]-self.variables.storage_charge[i,t] for i in heat_storage) + gb.quicksum(self.variables.Q[i,t] for i in CHP+heat_only+heat_pump),
                    gb.GRB.EQUAL,
                    self.data.heat_load[t],name='heat balance({0})'.format(t))    
                  
        # storage (1st stage)

        # heat storage: storage states update

        self.constraints.storage_update={}
        self.constraints.storage_init={}
        self.constraints.storage_final={}
        
        for h in heat_storage:
            
            for (t1,t2) in zip(time[:-1],time[1:]):
                self.constraints.storage_update[h,t2]=m.addConstr(
                    self.variables.storage_energy[h,t2],
                    gb.GRB.EQUAL,
                    self.variables.storage_energy[h,t1]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,t2]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,t2]-self.data.storage_loss[h])
        

            self.constraints.storage_init[h]=m.addConstr(
                self.variables.storage_energy[h,time[0]],
                gb.GRB.EQUAL,
                self.data.storage_energy_init[h]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,time[0]]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,time[0]]-self.data.storage_loss[h])


            self.constraints.storage_final[h]=m.addConstr(
                self.variables.storage_energy[h,time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.data.storage_energy_init[h])
                  
#%%

Q_sequential = {}
elec_maxprod_sequential = {}
elec_minprod_sequential = {}
        
for W_max in W_range:

    heat_dispatch = sequential_heat(W_max)                  
    heat_dispatch.model.params.OutputFlag = 0
    heat_dispatch.optimize()
    
    for t in time:
        
        for g in CHP+heat_only+heat_pump:
            Q_sequential[W_max,g,t]=heat_dispatch.variables.Q[g,t].x

        for g in heat_storage:
            Q_sequential[W_max,g,t]=heat_dispatch.variables.storage_discharge[g,t].x - heat_dispatch.variables.storage_charge[g,t].x   

    for t in time:
        
        for h in CHP_sorted['bp']:
            elec_maxprod_sequential[W_max,h,t] = 0
    
        for h in CHP_sorted['ex']:
            elec_maxprod_sequential[W_max,h,t] = (CHP_maxprod[h]-rho_heat[h]*Q_sequential[W_max,h,t])/rho_elec[h] - r_min[h]*Q_sequential[W_max,h,t]    

        for g in heat_pump:
            elec_maxprod_sequential[W_max,g,t] = Q_sequential[W_max,g,t]/COP[g]
            
        for g in gen:
            elec_maxprod_sequential[W_max,g,t] = elec_maxprod[g]   

        for g in wind:
            elec_maxprod_sequential[W_max,g,t] = W_max 
            
        for h in CHP:
            elec_minprod_sequential[W_max,h,t] = r_min[h]*Q_sequential[W_max,h,t]
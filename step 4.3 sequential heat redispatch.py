#%%
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:43:46 2016

@author: lemitri
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:31:11 2016

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
         

#%% building the STOCHSTICA MPEC optimization problem (GUROBI) = heat market clearing

class expando(object):
    '''
    
    
        A small class which can have attributes set
    '''
    pass

class sequential_heat_redispatch_scenario:
    def __init__(self,W_max,s0):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data(W_max,s0)
        self._build_model()

    
    def optimize(self):
        self.model.optimize()
    
    def _load_data(self,W_max,s0):



        #indexes
        self.data.time = time
        self.data.time_list=time_list
        self.data.gen=gen
        self.data.heat_storage=heat_storage
        self.data.heat_pump =heat_pump
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP
        
        # LOADS
        self.data.heat_load = heat_load
        
        # Heat station parameters
        self.data.CHP_maxprod = CHP_maxprod
        self.data.heat_maxprod = heat_maxprod
        self.data.rho_elec = rho_elec
        self.data.rho_heat = rho_heat
        self.data.r_min = r_min
        self.data.storage_discharge_eff= storage_rho_plus
        self.data.storage_charge_eff= storage_rho_moins
        self.data.storage_maxcapacity= storage_maxcapacity
        self.data.storage_maxprod= storage_maxprod
        self.data.storage_loss= storage_loss
        self.data.storage_energy_init= storage_init
        self.data.COP = COP 
        
        # Heat and electricity initial dispatch (DA)
        self.data.Q = {(h,t): Q_sequential[W_max,h,t] for h in CHP+heat_only for t in time}
        self.data.P = {(h,t):P_sequential_scenario[W_max,s0,h,t] for h in CHP+heat_pump for t in time}

        # Cost parameters

        self.data.alpha_up = {(g,t):alpha_up_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP+heat_only for t in time}       
        self.data.alpha_down = {(g,t):alpha_down_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP+heat_only for t in time}       

   
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.data.time
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP
        m = self.model
 
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

        self.variables.Q_up = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_up[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h],name='Q down({0},{1})'.format(h,t))

        self.variables.Q_down = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_down[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h],name='Q down({0},{1})'.format(h,t))

        m.update()
    
    def _build_objective(self): # building the objective function for the heat maret clearing

        #indexes shortcuts 
        time = self.data.time
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP
        m = self.model      

        m.setObjective(gb.quicksum(self.data.alpha_up[h,t]*self.variables.Q_up[h,t] -self.data.alpha_down[h,t]*self.variables.Q_down[h,t] for h in CHP+heat_only for t in time),   
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

        # heat balance
 
        self.constraints.heat_balance = {} 
        
        for t in time:
            self.constraints.heat_balance[t] = m.addConstr(
                    gb.quicksum(self.data.Q[h,t] + self.variables.Q_up[h,t] - self.variables.Q_down[h,t] for h in CHP+heat_only)+gb.quicksum(self.data.P[h,t]*self.data.COP[h] for h in heat_pump)+gb.quicksum(self.variables.storage_discharge[h,t]-self.variables.storage_charge[h,t] for h in heat_storage),
                    gb.GRB.EQUAL,
                    self.data.heat_load[t],name='heat balance ({0})'.format(t))


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


        #CHP's joint FOR in each scenario 
        
        self.constraints.heat_minprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_minprod[h,t] = m.addConstr(
                    self.data.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

        self.constraints.heat_maxprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_maxprod[h,t] = m.addConstr(
                    self.data.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t],
                    gb.GRB.LESS_EQUAL,
                    self.data.heat_maxprod[h])
                        
        self.constraints.CHP_maxprod = {} 
        self.constraints.CHP_ratio = {}
        
        for t in time:
            for h in CHP_sorted['ex']:
                
                self.constraints.CHP_maxprod[h,t] = m.addConstr(
                    self.data.rho_heat[h]*(self.data.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t])+self.data.rho_elec[h]*self.data.P[h,t],
                    gb.GRB.LESS_EQUAL,
                    self.data.CHP_maxprod[h])
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    self.data.P[h,t],
                    gb.GRB.GREATER_EQUAL,
                    self.data.r_min[h]*(self.data.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t]))

            for h in CHP_sorted['bp']:
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    self.data.P[h,t],
                    gb.GRB.EQUAL,
                    self.data.r_min[h]*(self.data.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t]))
                   

#%%

syst_cost_sequential_scenario = {}
heat_cost_sequential_scenario = {}
heat_price_sequential_scenario = {}

heat_redispatch_cost_sequential_average={}
heat_redispatch_cost_sequential_scenario={}

Q_sequential_scenario = {}
Q_up_sequential_scenario = {}
Q_down_sequential_scenario = {}
E_sequential_scenario = {}

Q_sequential_average = {}

syst_cost_sequential_average = {}
heat_cost_sequential_average = {}
heat_price_sequential_average = {}

for W_max in W_range:   
    
    for s0 in scenario_kmeans:
        dispatch = sequential_heat_redispatch_scenario(W_max,s0)
        dispatch.model.params.OutputFlag = 0
        dispatch.optimize()
           
        for g in heat_only+CHP:
            for t in time:
                Q_sequential_scenario[W_max,s0,g,t] = Q_sequential[W_max,g,t] + dispatch.variables.Q_up[g,t].x -  dispatch.variables.Q_down[g,t].x
                Q_up_sequential_scenario[W_max,s0,g,t] = dispatch.variables.Q_up[g,t].x 
                Q_down_sequential_scenario[W_max,s0,g,t] = dispatch.variables.Q_down[g,t].x
                
        for g in heat_pump:
            for t in time:
                Q_sequential_scenario[W_max,s0,g,t] = P_sequential_scenario[W_max,s0,g,t]*COP[g]
                
        for g in heat_storage:
            for t in time:
                Q_sequential_scenario[W_max,s0,g,t]=dispatch.variables.storage_discharge[g,t].x - dispatch.variables.storage_charge[g,t].x     
                
            E_sequential_scenario[W_max,s0,g,time[0]] = storage_init[g]
            
            for (t1,t2) in zip(time[:-1],time[1:]):
                E_sequential_scenario[W_max,s0,g,t2]=  dispatch.variables.storage_energy[g,t1].x  
    
        for t in time:
            heat_price_sequential_scenario[W_max,s0,t] = dispatch.constraints.heat_balance[t].Pi



for W_max in W_range:   
    
    for s0 in scenario_kmeans:
         
        heat_cost_sequential_scenario[W_max,s0] = sum(Q_sequential[W_max,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in heat_only for t in time) + sum((P_sequential_scenario[W_max,s0,g,t]*rho_elec[g]+Q_sequential[W_max,g,t]*rho_heat[g])*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] - P_sequential_scenario[W_max,s0,g,t]*spot_price_sequential_scenario[W_max,s0,t] for g in CHP for t in time) + sum(Q_sequential_scenario[W_max,s0,g,t]/COP[g]*spot_price_sequential_scenario[W_max,s0,t] for g in heat_pump for t in time)
        heat_redispatch_cost_sequential_scenario[W_max,s0] = sum(Q_up_sequential_scenario[W_max,s0,g,t]*alpha_up_kmeans[scenario_dic_kmeans[s0][2],g,t] - Q_down_sequential_scenario[W_max,s0,g,t]*alpha_down_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP+heat_only for t in time)
        syst_cost_sequential_scenario[W_max,s0] = sum(P_sequential_scenario[W_max,s0,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in gen for t in time) + sum(Q_sequential_scenario[W_max,s0,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in heat_only for t in time) + sum((P_sequential_scenario[W_max,s0,g,t]*rho_elec[g]+Q_sequential_scenario[W_max,s0,g,t]*rho_heat[g])*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP for t in time)
       
    for t in time:
        heat_price_sequential_average[W_max,t] = sum(heat_price_sequential_scenario[W_max,s,t] for s in scenario_kmeans)/S_all_kmeans
        for h in CHP+heat_only+heat_pump+heat_storage:
            Q_sequential_average[W_max,h,t] = sum(Q_sequential_scenario[W_max,s,h,t] for s in scenario_kmeans)/S_all_kmeans
    syst_cost_sequential_average[W_max] = sum(syst_cost_sequential_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    heat_cost_sequential_average[W_max] = sum(heat_cost_sequential_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    heat_redispatch_cost_sequential_average[W_max] = sum(heat_redispatch_cost_sequential_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans    
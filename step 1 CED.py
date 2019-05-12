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

W_range=[50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500]

    
#%% building the deterministic combined economic dispatch (CED) for each realization of the uncertain parameters (wind, loads, supply curve)
# each optimization problem is deterministic because the sources of uncertainty are revealed
    
class expando(object):
    '''
    
    
        A small class which can have attributes set
    '''
    pass

class CED_scenario:
    def __init__(self,s0,W_max):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data(s0,W_max)
        self._build_model()

    
    def optimize(self):
        self.model.optimize()

    def computeIIS(self):
        self.model.computeIIS()
    
    def _load_data(self,s0,W_max):
        
        #indexes
        self.data.time = time
        self.data.time_list=time_list
        self.data.gen=gen
        self.data.heat_storage=heat_storage
        self.data.elec_storage=elec_storage
        self.data.heat_pump =heat_pump
        self.data.wind=wind
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
        
        self.data.elec_maxprod = elec_maxprod
        
        self.data.W_max = W_max
        
        self.data.elec_load = {t:elec_load_scenario_kmeans[scenario_dic_kmeans[s0][0],t] for t in time}  
        self.data.wind_scenario = {(w,t):wind_scenario_kmeans[scenario_dic_kmeans[s0][1],w,t] for w in wind for t in time} #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!
        self.data.alpha = {(g,t):alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in gen+CHP+heat_only+wind for t in time} #BIDS

   
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.data.time
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP
      
        m = self.model
        
        #heat market optimization variables

        self.variables.Q = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h],name='Q({0},{1})'.format(h,t))

        self.variables.storage_discharge = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage+elec_storage:
                self.variables.storage_discharge[h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage plus({0},{1})'.format(h,t))
                    
        self.variables.storage_charge = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage+elec_storage:
                self.variables.storage_charge[h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage moins({0},{1})'.format(h,t))
                    

        self.variables.storage_energy = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage+elec_storage:
                self.variables.storage_energy[h,t] = m.addVar(lb=0,ub=self.data.storage_maxcapacity[h],name='storage energy({0},{1})'.format(h,t))                

        #electricity market optimization variables : primal variables
        
        self.variables.P = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP+gen+wind:
                self.variables.P[g,t] = m.addVar(lb=0,ub=self.data.elec_maxprod[g],name='P({0},{1})'.format(g,t)) # dispatch of electricity generators
            
        self.variables.HP_max_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_max_load[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h]/self.data.COP[h],name='HP max load (desired heat production)({0},{1})'.format(h,t))                

        m.update()
    
    def _build_objective(self): # building the objective function for the heat maret clearing

        #indexes shortcuts 
        time = self.data.time
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP      
        m = self.model     
              
        m.setObjective(gb.quicksum(self.data.alpha[g,t]*self.variables.P[g,t] for t in time for g in gen+wind)+gb.quicksum(self.data.alpha[g,t]*self.variables.Q[g,t] for t in time for g in heat_only)+gb.quicksum(self.data.alpha[g,t]*(self.data.rho_elec[g]*self.variables.P[g,t]+self.data.rho_heat[g]*self.variables.Q[g,t]) for t in time for g in CHP),   
            gb.GRB.MINIMIZE)
            
        
    def _build_constraints(self):

        #indexes shortcuts 
        time = self.data.time
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP     
        m = self.model
        
        self.constraints.heat_balance = {}
        
        for t in time:
                
                self.constraints.heat_balance[t] = m.addConstr(
                    gb.quicksum(self.variables.Q[i,t] for i in CHP+heat_only)+ gb.quicksum(self.data.COP[i]*self.variables.HP_max_load[i,t] for i in heat_pump)+gb.quicksum(self.variables.storage_discharge[h,t]-self.variables.storage_charge[h,t] for h in heat_storage),
                    gb.GRB.EQUAL,
                    self.data.heat_load[t],name='heat balance({0})'.format(t))    

        ##1) CHPs
                    
        self.constraints.CHP_maxprod = {} 
        self.constraints.CHP_ratio = {}
        
        for t in time:
                for h in self.data.CHP_sorted['ex']:
                    
                    self.constraints.CHP_maxprod[h,t] = m.addConstr(
                        self.data.rho_heat[h]*self.variables.Q[h,t]+self.data.rho_elec[h]*self.variables.P[h,t],
                        gb.GRB.LESS_EQUAL,
                        self.data.CHP_maxprod[h],name='CHP maxprod({0},{1})'.format(h,t))
                    
                    self.constraints.CHP_ratio[h,t] = m.addConstr(
                        self.variables.P[h,t],
                        gb.GRB.GREATER_EQUAL,
                        self.data.r_min[h]*self.variables.Q[h,t],name='CHP ratio({0},{1})'.format(h,t))

                for h in self.data.CHP_sorted['bp']:
                    
                    self.constraints.CHP_ratio[h,t] = m.addConstr(
                        self.variables.P[h,t],
                        gb.GRB.EQUAL,
                        self.data.r_min[h]*self.variables.Q[h,t],name='CHP ratio({0},{1})'.format(h,t)) 
        
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
                            

        # wind realization

        self.constraints.wind_scenario = {}
        
        for t in time:
                for g in wind: 
                    
                    self.constraints.wind_scenario[g,t] = m.addConstr(
                        self.variables.P[g,t],
                        gb.GRB.LESS_EQUAL,
                        self.data.wind_scenario[g,t]*self.data.W_max,name='wind scenario({0},{1})'.format(g,t))
        

        self.constraints.elec_balance = {}
        
        for t in time:
                    
                self.constraints.elec_balance[t] = m.addConstr(
                    gb.quicksum(self.variables.P[g,t] for g in gen+wind+CHP),
                    gb.GRB.EQUAL,
                    self.data.elec_load[t]+gb.quicksum(self.variables.HP_max_load[g,t] for g in heat_pump),name='elec balance({0})'.format(t))   


#%% SOLVE FOR EACH REALIZATION OF SCENARIOS and different wind power penetration (W_max)

syst_cost_CED_scenario = {}
heat_cost_CED_scenario = {}
elec_cost_CED_scenario = {}
heat_price_CED_scenario = {}
spot_price_CED_scenario = {}

P_CED_scenario  = {}
Q_CED_scenario = {}
E_CED_scenario = {}

wind_curtailment_CED_scenario = {}

Q_CED_average = {}

heat_cost_CED_average = {}
elec_cost_CED_average = {}
heat_price_CED_average = {}
spot_price_CED_average = {}
syst_cost_CED_average = {}
wind_curtailment_CED_average = {}

for W_max in W_range:
    for s0 in scenario_kmeans:
        dispatch_CED = CED_scenario(s0,W_max)
        dispatch_CED.model.params.OutputFlag = 0
        dispatch_CED.optimize()
    
        for g in gen+wind+CHP:
            for t in time:
                P_CED_scenario[W_max,s0,g,t]=dispatch_CED.variables.P[g,t].x
    
           
        for g in heat_only+CHP:
            for t in time:
                Q_CED_scenario[W_max,s0,g,t]=dispatch_CED.variables.Q[g,t].x
        for g in heat_pump:
            for t in time:
                Q_CED_scenario[W_max,s0,g,t]=COP[g]*dispatch_CED.variables.HP_max_load[g,t].x
        for g in heat_storage:
            for t in time:
                Q_CED_scenario[W_max,s0,g,t]=dispatch_CED.variables.storage_discharge[g,t].x - dispatch_CED.variables.storage_charge[g,t].x     
            E_CED_scenario[W_max,s0,g,time[0]] = storage_init[g]
            for (t1,t2) in zip(time[:-1],time[1:]):
                E_CED_scenario[W_max,s0,g,t2]=  dispatch_CED.variables.storage_energy[g,t1].x  
            
        for t in time:
            spot_price_CED_scenario[W_max,s0,t] = dispatch_CED.constraints.elec_balance[t].Pi
            heat_price_CED_scenario[W_max,s0,t] = dispatch_CED.constraints.heat_balance[t].Pi

        wind_curtailment_CED_scenario[W_max,s0] = sum(wind_scenario_kmeans[scenario_dic_kmeans[s0][1],g,t]*W_max-P_CED_scenario[W_max,s0,g,t] for g in wind for t in time)
        syst_cost_CED_scenario[W_max,s0] = dispatch_CED.model.objval

        
for W_max in W_range:
    for s0 in scenario_kmeans:
        elec_cost_CED_scenario[W_max,s0] = sum(P_CED_scenario[W_max,s0,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in gen for t in time) + sum(P_CED_scenario[W_max,s0,g,t]*rho_elec[g]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP for t in time) 
        heat_cost_CED_scenario[W_max,s0] = sum(Q_CED_scenario[W_max,s0,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in heat_only for t in time) + sum((P_CED_scenario[W_max,s0,g,t]*rho_elec[g] + Q_CED_scenario[W_max,s0,g,t]*rho_heat[g])*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] - P_CED_scenario[W_max,s0,g,t]*spot_price_CED_scenario[W_max,s0,t] for g in CHP for t in time) + sum(Q_CED_scenario[W_max,s0,g,t]/COP[g]*spot_price_CED_scenario[W_max,s0,t] for g in heat_pump for t in time)
        
for W_max in W_range:  
    for t in time:
        spot_price_CED_average[W_max,t] = sum(spot_price_CED_scenario[W_max,s,t] for s in scenario_kmeans)/S_all_kmeans
        heat_price_CED_average[W_max,t] = sum(heat_price_CED_scenario[W_max,s,t] for s in scenario_kmeans)/S_all_kmeans

        for h in CHP+heat_only+heat_pump+heat_storage:
            Q_CED_average[W_max,h,t]=sum(Q_CED_scenario[W_max,s0,g,t] for s0 in scenario_kmeans)/S_all_kmeans
            
    syst_cost_CED_average[W_max] = sum(syst_cost_CED_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    wind_curtailment_CED_average[W_max] = sum(wind_curtailment_CED_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    elec_cost_CED_average[W_max] = sum(elec_cost_CED_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    heat_cost_CED_average[W_max] = sum(heat_cost_CED_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans

for W_max in W_range:  
    for t in time:
        for h in CHP+heat_only+heat_pump+heat_storage:
            Q_CED_average[W_max,h,t]=sum(Q_CED_scenario[W_max,s0,h,t] for s0 in scenario_kmeans)/S_all_kmeans
    
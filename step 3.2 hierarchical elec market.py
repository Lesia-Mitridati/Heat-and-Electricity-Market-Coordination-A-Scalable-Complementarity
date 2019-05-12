# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:51:31 2016

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

#%% NORDPOOL MC clearing for each scenario IN SAMPLE !!!!

class expando(object):
    '''
        A small class which can have attributes set
    '''
    pass

class hierarchical_elec_scenario:
    def __init__(self,s0,W_max):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data(s0,W_max)
        self._build_model()
    
    def optimize(self):
        self.model.optimize()
    
    def _load_data(self,s0,W_max):
        
        
        #indexes
        self.data.W_max = W_max
        self.data.time = time
        self.data.time_list=time_list
        self.data.gen=gen
        self.data.heat_pump =heat_pump
        self.data.wind=wind
        self.data.CHP = CHP
        self.data.rho_elec = rho_elec
        
        # scenarios
        self.data.elec_load = {t:elec_load_scenario_kmeans[scenario_dic_kmeans[s0][0],t] for t in time}  #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!            
        self.data.wind_scenario = {(w,t):wind_scenario_kmeans[scenario_dic_kmeans[s0][1],w,t] for w in wind for t in time} #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!
        self.data.alpha = {(g,t):alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP+gen+wind for t in time} #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!
        
        # Elec station parameters
        self.data.elec_maxprod = {(g,t):elec_maxprod_hierarchical[W_max,g,t] for g in gen+wind+CHP+heat_pump for t in time}
        self.data.elec_minprod = {(g,t):elec_minprod_hierarchical[W_max,g,t] for g in CHP for t in time}

        self.data.alpha_HP = {(g,t):100 for g in heat_pump for t in time}

    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
    
    def _build_variables(self):
        
        #indexes
        m = self.model
        time = self.data.time
        time_list = self.data.time_list
        gen = self.data.gen
        heat_pump = self.data.heat_pump 
        wind = self.data.wind
        CHP = self.data.CHP 

        #electricity market optimization variables : primal variables
        
        self.variables.P = {} # electricity production from electricity generators
        for t in time:
            for g in CHP+gen+wind:
                self.variables.P[g,t] = m.addVar(lb=0,ub=self.data.elec_maxprod[g,t],name='elec prod at marginal cost({0},{1})'.format(g,t)) # dispatch of electricity generators


        self.variables.HP_load = {} # electricity production from electricity generators
        for t in time:
            for g in heat_pump:
                self.variables.HP_load[g,t] = m.addVar(lb=0,ub=self.data.elec_maxprod[g,t],name='elec prod at marginal cost({0},{1})'.format(g,t))
                
                
        self.variables.P_min = {} # electricity production from electricity generators
        for t in time:
            for g in CHP:
                self.variables.P_min[g,t] = m.addVar(lb=0,ub=self.data.elec_minprod[g,t],name='elec MIN prod for CHPs({0},{1})'.format(g,t)) # dispatch of electricity generators
  
        m.update()
        
        
    def _build_objective(self): # building the objective function for the heat maret clearing
        
        #indexes
        m = self.model
        time = self.data.time
        time_list = self.data.time_list
        gen = self.data.gen
        wind = self.data.wind
        CHP = self.data.CHP 
        
        m.setObjective(-gb.quicksum(self.data.alpha_HP[g,t]*self.variables.HP_load[g,t] for t in time for g in heat_pump)+gb.quicksum(self.data.alpha[g,t]*self.variables.P[g,t] for t in time for g in gen+wind)+gb.quicksum(self.data.alpha[g,t]*self.data.rho_elec[g]*self.variables.P[g,t]-300*self.variables.P_min[g,t] for t in time for g in CHP),
            gb.GRB.MINIMIZE)
            
        
    def _build_constraints(self):
        
        #indexes
        m = self.model
        time = self.data.time
        time_list = self.data.time_list
        gen = self.data.gen
        heat_pump = self.data.heat_pump 
        wind = self.data.wind
        CHP = self.data.CHP 
 
        # wind realization

        self.constraints.wind_scenario = {}
        
        for t in time:
                for g in wind: 
                    
                    self.constraints.wind_scenario[g,t] = m.addConstr(
                        self.variables.P[g,t],
                        gb.GRB.LESS_EQUAL,
                        self.data.wind_scenario[g,t]*self.data.elec_maxprod[g,t],name='wind scenario({0},{1})'.format(g,t))
        

        self.constraints.elec_balance = {}
        
        for t in time:
                    
                self.constraints.elec_balance[t] = m.addConstr(
                    gb.quicksum(self.variables.P_min[g,t] for g in CHP)+gb.quicksum(self.variables.P[g,t] for g in CHP+gen+wind),
                    gb.GRB.EQUAL,
                    self.data.elec_load[t]+gb.quicksum(self.variables.HP_load[g,t] for g in heat_pump),name='elec balance({0})'.format(t)) 


#%% solve the MC for different realisations of the scenarios 

elec_cost_hierarchical_scenario = {}
spot_price_hierarchical_scenario = {}

P_hierarchical_scenario  = {}

wind_curtailment_hierarchical_scenario = {}

elec_cost_hierarchical_average = {}
spot_price_hierarchical_average = {}
P_hierarchical_average = {}
wind_curtailment_hierarchical_average = {}

for W_max in W_range:   
    
    for s0 in scenario_kmeans:
        
        dispatch = hierarchical_elec_scenario(s0,W_max)
        dispatch.model.params.OutputFlag = 0
        dispatch.optimize()

    
        for g in gen+wind:
            for t in time:
                
                P_hierarchical_scenario[W_max,s0,g,t]=dispatch.variables.P[g,t].x
    
        for g in CHP:
            for t in time:
                
                P_hierarchical_scenario[W_max,s0,g,t]=dispatch.variables.P_min[g,t].x + dispatch.variables.P[g,t].x
        
        for g in heat_pump:
            for t in time:
                P_hierarchical_scenario[W_max,s0,g,t]=dispatch.variables.HP_load[g,t].x
                
        
        for t in time:

            spot_price_hierarchical_scenario[W_max,s0,t] = dispatch.constraints.elec_balance[t].Pi

        wind_curtailment_hierarchical_scenario[W_max,s0] = sum(W_max*wind_scenario_kmeans[scenario_dic_kmeans[s0][1],w,t]- P_hierarchical_scenario[W_max,s0,w,t] for w in wind for t in time)

for W_max in W_range:   
    
    for s0 in scenario_kmeans:
        elec_cost_hierarchical_scenario[W_max,s0] = sum(P_hierarchical_scenario[W_max,s0,g,t]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in gen for t in time) + sum(P_hierarchical_scenario[W_max,s0,g,t]*rho_elec[g]*alpha_scenario_kmeans[scenario_dic_kmeans[s0][2],g,t] for g in CHP for t in time) 

    for t in time:
        for h in CHP+gen+wind+heat_pump:
            
            P_hierarchical_average[W_max,g,t]=sum(P_hierarchical_scenario[W_max,s0,g,t] for s0 in scenario_kmeans)/S_all_kmeans
            
        
        spot_price_hierarchical_average[W_max,t] = sum(spot_price_hierarchical_scenario[W_max,s,t] for s in scenario_kmeans)/S_all_kmeans

    wind_curtailment_hierarchical_average[W_max] = sum(wind_curtailment_hierarchical_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans
    elec_cost_hierarchical_average[W_max] = sum(elec_cost_hierarchical_scenario[W_max,s] for s in scenario_kmeans)/S_all_kmeans

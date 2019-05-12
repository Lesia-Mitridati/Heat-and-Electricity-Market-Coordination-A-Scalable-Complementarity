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
import time as ttt
            
class expando(object):
    '''
    
    
        A small class which can have attributes set
    '''
    pass

#%%

class hierarchical_heat:
    def __init__(self,W_max):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data(W_max)
        self._build_model()

    
    def optimize(self):
        self.model.optimize()


    def computeIIS(self):
        self.model.computeIIS()
        
        
    def _load_data(self,W_max):

        #indexes
        self.data.S_kmeans_0 = S_kmeans_0
        self.data.S_all_kmeans_0 = S_all_kmeans_0
        self.data.scenario_load_kmeans_0=scenario_load_kmeans_0
        self.data.scenario_wind_kmeans_0=scenario_wind_kmeans_0
        self.data.scenario_supply_kmeans_0=scenario_supply_kmeans_0
        self.data.scenario_kmeans_0=scenario_kmeans_0
        self.data.scenario_list_kmeans_0=scenario_list_kmeans_0
        self.data.scenario_dic_kmeans_0=scenario_dic_kmeans_0
    
        self.data.time = time
        self.data.gen=gen
        self.data.heat_storage=heat_storage
        self.data.elec_storage=elec_storage
        self.data.heat_pump =heat_pump
        self.data.wind=wind
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP  
        
        # Cost parameters
        self.data.alpha_min = -300
        self.data.penalty  = 0.001 #PENALTY for wind curtailment  
        self.data.alpha_HP = 100

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
#        self.data.alpha_HP = {h:max([alpha_m[h]*rho_elec[h] for h in CHP])*COP[h] for h in heat_pump}
        # LOADS
        self.data.heat_load = heat_load
        
        # Elec station parameters
        self.data.elec_maxprod = elec_maxprod
        self.data.W_max = W_max

        self.data.elec_load = {(s,t):elec_load_scenario_kmeans_0[self.data.scenario_dic_kmeans_0[s][0],t] for s in self.data.scenario_kmeans_0 for t in time}  
        self.data.wind_scenario = {(s,w,t):wind_scenario_kmeans_0[self.data.scenario_dic_kmeans_0[s][1],w,t] for s in self.data.scenario_kmeans_0 for w in wind for t in time} #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!
        self.data.alpha = {(s,g,t):alpha_scenario_kmeans_0[self.data.scenario_dic_kmeans_0[s][2],g,t] for g in CHP+heat_only+gen+wind for s in self.data.scenario_kmeans_0 for t in time} #BIDS
        self.data.alpha_up = {(s,g,t):alpha_up_kmeans_0[self.data.scenario_dic_kmeans_0[s][2],g,t] for g in CHP+heat_only for s in self.data.scenario_kmeans_0 for t in time}
        self.data.alpha_down = {(s,g,t):alpha_down_kmeans_0[self.data.scenario_dic_kmeans_0[s][2],g,t] for g in CHP+heat_only for s in self.data.scenario_kmeans_0 for t in time}       
   
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.data.time
        scenario_kmeans_0 = self.data.scenario_kmeans_0
        S_all_kmeans_0 = self.data.S_all_kmeans_0
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        #elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP     
        m = self.model
        
        # heat market optimization variables

        ## 1st stage variables 

        self.variables.Q = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q[h,t] = m.addVar(lb=0.00,ub=self.data.heat_maxprod[h],name='Q({0},{1})'.format(h,t))


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

        self.variables.HP_max_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_max_load[h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h]/self.data.COP[h],name='HP max load (desired heat production)({0},{1})'.format(h,t))                
                
                
        ## 2nd stage variables

        self.variables.Q_down = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                for s in scenario_kmeans_0:
                    self.variables.Q_down[s,h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h],name='Q down({0},{1},{2})'.format(s,h,t))

        self.variables.Q_up = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                for s in scenario_kmeans_0:
                    self.variables.Q_up[s,h,t] = m.addVar(lb=0,ub=self.data.heat_maxprod[h],name='Q up({0},{1},{2})'.format(s,h,t))
           
        self.variables.storage_discharge_scenario = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                for s in scenario_kmeans_0:
                    self.variables.storage_discharge_scenario[s,h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage discharge scenario({0},{1},{2})'.format(s,h,t))
                    
        self.variables.storage_charge_scenario = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                for s in scenario_kmeans_0:
                    self.variables.storage_charge_scenario[s,h,t] = m.addVar(lb=0,ub=self.data.storage_maxprod[h],name='storage charge scenario({0},{1},{2})'.format(s,h,t))
                    

        self.variables.storage_energy_scenario = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                for s in scenario_kmeans_0:
                    self.variables.storage_energy_scenario[s,h,t] = m.addVar(lb=0,ub=self.data.storage_maxcapacity[h],name='storage energy scenario({0},{1},{2})'.format(s,h,t))                

        # electricity market optimization variables

        ## primal variables

        self.variables.P = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in gen+wind+CHP_sorted['ex']:
                for s in scenario_kmeans_0:
                    self.variables.P[s,g,t] = m.addVar(lb=0,name='P({0},{1},{2})'.format(s,g,t)) # dispatch of electricity generators


        self.variables.P_0 = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP:
                for s in scenario_kmeans_0:
                    self.variables.P_0[s,g,t] = m.addVar(lb=0,name='P 0({0},{1},{2})'.format(s,g,t)) # dispatch of electricity generators              


        self.variables.HP_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                for s in scenario_kmeans_0:
                    self.variables.HP_load[s,h,t] = m.addVar(lb=0,name='HP max load (desired heat production)({0},{1})'.format(h,t))   
         
                    
        ## dual variables
        
        self.variables.lambda_powerbalance = {}
        for t in time:
            for s in scenario_kmeans_0:
                self.variables.lambda_powerbalance[s,t] = m.addVar(lb=-gb.GRB.INFINITY,name='lambda powerbalance({0},{1})'.format(s,t))
        
        self.variables.mu_max = {}
        for t in time:
            for g in CHP_sorted['ex']+gen+wind+heat_pump:
                for s in scenario_kmeans_0:
                    self.variables.mu_max[s,g,t] = m.addVar(lb=0,name='mu max({0},{1},{2})'.format(s,g,t))
                    
        self.variables.mu_min = {}
        for t in time:
            for g in CHP_sorted['ex']+gen+wind+heat_pump:
                for s in scenario_kmeans_0:
                    self.variables.mu_min[s,g,t] = m.addVar(lb=0,name='mu min({0},{1},{2})'.format(s,g,t))


        self.variables.mu_max_0 = {}
        for t in time:
            for g in CHP:
                for s in scenario_kmeans_0:
                    self.variables.mu_max_0[s,g,t] = m.addVar(lb=0,name='mu max 0({0},{1},{2})'.format(s,g,t))
                    
        self.variables.mu_min_0 = {}
        for t in time:
            for g in CHP:
                for s in scenario_kmeans_0:
                    self.variables.mu_min_0[s,g,t] = m.addVar(lb=0,name='mu min 0({0},{1},{2})'.format(s,g,t))
                    
                    
                    
        self.variables.f_scenario = {}
        for s in scenario_kmeans_0:
            self.variables.f_scenario[s] = m.addVar(lb=-gb.GRB.INFINITY,name='obj value for scenario({0})'.format(s))


        ## SOS1 variables (POSITIVE)


        self.variables.SOS_max = {} 
        
        for t in time:
            for h in CHP_sorted['ex']+gen+wind+heat_pump:
                for s in scenario_kmeans_0:
                    self.variables.SOS_max[s,h,t] = m.addVar(lb=0,name='SOS max({0},{1},{2})'.format(s,h,t))        
        
        self.variables.SOS_min = {} 
        
        for t in time:
            for h in CHP_sorted['ex']+gen+wind+heat_pump:
                for s in scenario_kmeans_0:
                    self.variables.SOS_min[s,h,t] = m.addVar(lb=0,name='SOS min({0},{1},{2})'.format(s,h,t))               


        self.variables.SOS_max_0 = {} 
        
        for t in time:
            for h in CHP:
                for s in scenario_kmeans_0:
                    self.variables.SOS_max_0[s,h,t] = m.addVar(lb=0,name='SOS max 0({0},{1},{2})'.format(s,h,t))        
        
        self.variables.SOS_min_0 = {} 
        
        for t in time:
            for h in CHP:
                for s in scenario_kmeans_0:
                    self.variables.SOS_min_0[s,h,t] = m.addVar(lb=0,name='SOS min 0({0},{1},{2})'.format(s,h,t))
                    
                    

        m.update()
    
    def _build_objective(self): # building the objective function for the heat maret clearing

        #indexes shortcuts 
        time = self.data.time
        scenario_kmeans_0 = self.data.scenario_kmeans_0
        S_all_kmeans_0 = self.data.S_all_kmeans_0
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        #elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP      
        m = self.model      
              
#        m.setObjective(gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_heat[h]*self.variables.Q[h,t] for s in scenario for h in CHP+heat_only for t in time)/S + self.data.penalty*gb.quicksum(self.data.wind_scenario[s,w,t]*self.data.elec_maxprod[w]-self.variables.P[s,w,t] for w in wind for t in time for s in scenario)/S,   
#            gb.GRB.MINIMIZE)

        m.setObjective(gb.quicksum(self.variables.f_scenario[s] for s in scenario_kmeans_0)/S_all_kmeans_0,   
            gb.GRB.MINIMIZE)
            
#        m.setObjective(gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_heat[h]*self.variables.Q[h,t] + self.data.alpha_up[s,h,t]*self.data.rho_heat[h]*self.variables.Q_up[s,h,t] - self.data.alpha_down[s,h,t]*self.data.rho_heat[h]*self.variables.Q_down[s,h,t] for s in scenario for h in CHP+heat_only for t in time)/S_all + gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_elec[h]*self.variables.P[s,h,t] for s in scenario for h in CHP for t in time)/S_all + gb.quicksum(self.data.alpha[s,h,t]*self.variables.P[s,h,t] for s in scenario for h in gen+wind for t in time)/S_all + gb.quicksum(self.variables.mu_max[s,h,t]*self.data.elec_maxprod[h] for s in scenario for h in gen for t in time)/S_all + gb.quicksum(self.variables.mu_max[s,h,t]*self.data.wind_scenario[s,h,t]*self.data.elec_maxprod[h] for s in scenario for h in wind for t in time)/S_all - gb.quicksum(self.variables.lambda_powerbalance[s,t]*self.data.elec_load[s,t] for s in scenario for t in time)/S_all,   
#            gb.GRB.MINIMIZE)
        
    def _build_constraints(self):

        #indexes shortcuts 
        time = self.data.time
        scenario_kmeans_0 = self.data.scenario_kmeans_0
        S_all_kmeans_0 = self.data.S_all_kmeans_0
        gen=self.data.gen
        heat_storage=self.data.heat_storage
        #elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        wind=self.data.wind
        heat_only=self.data.heat_only
        CHP_sorted=self.data.CHP_sorted
        CHP=self.data.CHP      
        m = self.model

        # Upper Level Constraints
        # 1st stage constraints
        # heat balance

        self.constraints.obj_value = {} 
        
        for s in scenario_kmeans_0:
            self.constraints.obj_value[s] = m.addConstr(
                    self.variables.f_scenario[s],
                    gb.GRB.EQUAL,
                    gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_heat[h]*self.variables.Q[h,t] + self.data.alpha_up[s,h,t]*self.variables.Q_up[s,h,t] - self.data.alpha_down[s,h,t]*self.variables.Q_down[s,h,t] for h in CHP+heat_only for t in time) + gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_elec[h]*self.variables.P_0[s,h,t] for h in CHP for t in time) + gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_elec[h]*self.variables.P[s,h,t] for h in CHP_sorted['ex'] for t in time) + gb.quicksum(self.data.alpha[s,h,t]*self.variables.P[s,h,t] for h in gen+wind for t in time) + gb.quicksum(self.variables.mu_max[s,h,t]*self.data.elec_maxprod[h] for h in gen for t in time) + gb.quicksum(self.variables.mu_max[s,h,t]*self.data.wind_scenario[s,h,t]*self.data.W_max for h in wind for t in time) - gb.quicksum(self.variables.lambda_powerbalance[s,t]*self.data.elec_load[s,t] for t in time))

 
        self.constraints.heat_balance = {} 
        
        for t in time:
            self.constraints.heat_balance[t] = m.addConstr(
                    gb.quicksum(self.variables.Q[h,t] for h in CHP+heat_only)+gb.quicksum(self.data.COP[h]*self.variables.HP_max_load[h,t] for h in heat_pump)+gb.quicksum(self.variables.storage_discharge[h,t]-self.variables.storage_charge[h,t] for h in heat_storage),
                    gb.GRB.EQUAL,
                    self.data.heat_load[t],name='heat balance ({0})'.format(t))

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
                    self.variables.storage_energy[h,t1]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,t2]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,t2]-self.data.storage_loss[h],name='heat storage update({0},{1},{2})'.format(h,t1,t2))
        

            self.constraints.storage_init[h]=m.addConstr(
                self.variables.storage_energy[h,time[0]],
                gb.GRB.EQUAL,
                self.data.storage_energy_init[h]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,time[0]]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,time[0]]-self.data.storage_loss[h],name='heat storage init({0})'.format(h))


            self.constraints.storage_final[h]=m.addConstr(
                self.variables.storage_energy[h,time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.data.storage_energy_init[h],name='heat storage final({0})'.format(h))

        #2nd stage constraints

        #CHP's joint FOR in each scenario 
        
        self.constraints.heat_minprod = {} 

        for t in time:
            for s in scenario_kmeans_0:
                for h in CHP+heat_only:
                    
                    self.constraints.heat_minprod[s,h,t] = m.addConstr(
                        self.variables.Q[h,t]+self.variables.Q_up[s,h,t]-self.variables.Q_down[s,h,t],
                        gb.GRB.GREATER_EQUAL,
                        0,name='heat min prod({0},{1},{2})'.format(s,h,t))

        self.constraints.heat_maxprod = {} 

        for t in time:
            for s in scenario_kmeans_0:
                for h in CHP+heat_only:
                    
                    self.constraints.heat_maxprod[s,h,t] = m.addConstr(
                        self.variables.Q[h,t]+self.variables.Q_up[s,h,t]-self.variables.Q_down[s,h,t],
                        gb.GRB.LESS_EQUAL,
                        self.data.heat_maxprod[h],name='heat max prod({0},{1},{2})'.format(s,h,t))
                        
        self.constraints.CHP_maxprod = {} 
        self.constraints.CHP_ratio = {}
        
        for t in time:
            for s in scenario_kmeans_0:
                for h in CHP_sorted['ex']:
                    
                    self.constraints.CHP_maxprod[s,h,t] = m.addConstr(
                        self.data.rho_heat[h]*(self.variables.Q[h,t]+self.variables.Q_up[s,h,t]-self.variables.Q_down[s,h,t])+self.data.rho_elec[h]*(self.variables.P_0[s,h,t]+self.variables.P[s,h,t]),
                        gb.GRB.LESS_EQUAL,
                        self.data.CHP_maxprod[h],name='CHP max prod({0},{1},{2})'.format(s,h,t))
                    
                    self.constraints.CHP_ratio[s,h,t] = m.addConstr(
                        self.variables.P_0[s,h,t]+self.variables.P[s,h,t],
                        gb.GRB.GREATER_EQUAL,
                        self.data.r_min[h]*(self.variables.Q[h,t]+self.variables.Q_up[s,h,t]-self.variables.Q_down[s,h,t]),name='CHP ratio({0},{1},{2})'.format(s,h,t))

                for h in CHP_sorted['bp']:
                    
                    self.constraints.CHP_ratio[s,h,t] = m.addConstr(
                        self.variables.P_0[s,h,t],
                        gb.GRB.EQUAL,
                        self.data.r_min[h]*(self.variables.Q[h,t]+self.variables.Q_up[s,h,t]-self.variables.Q_down[s,h,t]),name='CHP ratio({0},{1},{2})'.format(s,h,t))


        self.constraints.heat_balance_scenario = {} 
        
        for t in time:
            for s in scenario_kmeans_0:
                self.constraints.heat_balance_scenario[s,t] = m.addConstr(
                        gb.quicksum(self.variables.Q[h,t] + self.variables.Q_up[s,h,t] - self.variables.Q_down[s,h,t] for h in CHP+heat_only)+gb.quicksum(self.data.COP[h]*self.variables.HP_load[s,h,t] for h in heat_pump)+gb.quicksum(self.variables.storage_discharge_scenario[s,h,t]-self.variables.storage_charge_scenario[s,h,t] for h in heat_storage),
                        gb.GRB.EQUAL,
                        self.data.heat_load[t],name='heat balance _scenario ({0},{1})'.format(s,t))
                    

        # storage (and stage)

        # heat storage: storage states update

        self.constraints.storage_update_scenario={}
        self.constraints.storage_init_scenario={}
        self.constraints.storage_final_scenario={}
        
        for s in scenario_kmeans_0:
            for h in heat_storage:        
                for (t1,t2) in zip(time[:-1],time[1:]):
                    self.constraints.storage_update_scenario[s,h,t2]=m.addConstr(
                        self.variables.storage_energy_scenario[s,h,t2],
                        gb.GRB.EQUAL,
                        self.variables.storage_energy_scenario[s,h,t1]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[s,h,t2]+self.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[s,h,t2]-self.data.storage_loss[h],name='storage update scenario ({0},{1},{2},{3})'.format(s,h,t1,t2))
    
                self.constraints.storage_init_scenario[h]=m.addConstr(
                    self.variables.storage_energy_scenario[s,h,time[0]],
                    gb.GRB.EQUAL,
                    self.data.storage_energy_init[h]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[s,h,time[0]]+self.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[s,h,time[0]]-self.data.storage_loss[h],name='storage init scenario({0},{1})'.format(s,h))
       
                self.constraints.storage_final_scenario[h]=m.addConstr(
                    self.variables.storage_energy_scenario[s,h,time[-1]],
                    gb.GRB.GREATER_EQUAL,
                    self.data.storage_energy_init[h],name='storage final scenario({0},{1})'.format(s,h))
                    
        # Lower level Constraints
        
        ## Equality
        
        self.constraints.elec_balance = {} 
        
        for t in time:
            for s in scenario_kmeans_0:
                self.constraints.elec_balance[s,t] = m.addConstr(
                    gb.quicksum(self.variables.P[s,g,t] for g in CHP_sorted['ex']+gen+wind)+gb.quicksum(self.variables.P_0[s,g,t] for g in CHP)-gb.quicksum(self.variables.HP_load[s,h,t] for h in heat_pump),
                    gb.GRB.EQUAL,
                    self.data.elec_load[s,t],name='elec balance({0},{1})'.format(s,t))       
                
#        ## Inequality
#        
#        self.constraints.P_max = {}
#        
#        for t in time:
#            for s in scenario:     
#                
#                for g in gen:
#                    self.constraints.P_max[s,g,t] = m.addConstr(
#                        self.variables.P[s,g,t],
#                        gb.GRB.LESS_EQUAL,
#                        self.data.elec_maxprod[g])
#           
#                for g in wind:                    
#                    self.constraints.P_max[s,g,t] = m.addConstr(
#                        self.variables.P[s,g,t],
#                        gb.GRB.LESS_EQUAL,
#                        self.data.wind_scenario[s,g,t]*self.data.elec_maxprod[g])
#
#                for g in CHP_sorted['ex']:
#                    self.constraints.P_max[s,g,t] = m.addConstr(
#                        self.variables.P[s,g,t],
#                        gb.GRB.LESS_EQUAL,
#                        (self.data.CHP_maxprod[g]-self.data.rho_heat[g]*self.variables.Q[g,t])/self.data.rho_elec[g])                           
#                                                
#                        
#                for g in CHP_sorted['bp']:
#                    self.constraints.P_max[s,g,t] = m.addConstr(
#                        self.variables.P[s,g,t],
#                        gb.GRB.LESS_EQUAL,
#                        self.data.r_min[g]*self.variables.Q[g,t])                                                       
         
        # lagrangian function: dL/dx + lambda.dh/dx + mu.dg/dx = 0
        
        self.constraints.dL = {} # gradient of lagragian function / dP_gen[g,t]
        
        for t in time:
            for s in scenario_kmeans_0:
                
                for g in gen+wind:
                    self.constraints.dL[s,g,t] = m.addConstr(
                        self.data.alpha[s,g,t]-self.variables.lambda_powerbalance[s,t]-self.variables.mu_min[s,g,t]+self.variables.mu_max[s,g,t],
                        gb.GRB.EQUAL,
                        0,name='lagrangian({0},{1},{2})'.format(s,g,t))
                                            
                for h in CHP_sorted['ex']:
                    self.constraints.dL['up',s,h,t] = m.addConstr(
                        self.data.alpha[s,h,t]*self.data.rho_elec[h]-self.variables.lambda_powerbalance[s,t]-self.variables.mu_min[s,h,t]+self.variables.mu_max[s,h,t],
                        gb.GRB.EQUAL,
                        0,name='lagangian({0},{1},{2})'.format(s,h,t))

                for h in CHP:
                    self.constraints.dL['0',s,h,t] = m.addConstr(
                        self.data.alpha_min-self.variables.lambda_powerbalance[s,t]-self.variables.mu_min_0[s,h,t]+self.variables.mu_max_0[s,h,t],
                        gb.GRB.EQUAL,
                        0,name='lagangian 0 ({0},{1},{2})'.format(s,h,t))


                for h in heat_pump:
                    self.constraints.dL[s,h,t] = m.addConstr(
                        -self.data.alpha_HP+self.variables.lambda_powerbalance[s,t]-self.variables.mu_min[s,h,t]+self.variables.mu_max[s,h,t],
                        gb.GRB.EQUAL,
                        0,name='lagangian 0 ({0},{1},{2})'.format(s,h,t))
                        
        # complementarity constraintes mu.g = 0 
        # reformulate as mu <= M*(1-u) and g <= M*u


        # SOS constraints 
        ## SOS.mu = 0
        ## SOS = g >= 0 


        self.constraints.SOS_min_1 = {}
        self.constraints.SOS_min_2 = {}

        for t in time:
            for g in gen+wind: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_min_1[s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_min[s,g,t],
                        self.variables.mu_min[s,g,t]])

                    self.constraints.SOS_min_2[s,g,t] = m.addConstr(
                        self.variables.SOS_min[s,g,t],
                        gb.GRB.EQUAL,
                        self.variables.P[s,g,t],
                        name='min elec production SOS 2 ({0},{1},{2})'.format(s,g,t))                        
                    
        for t in time:
            for g in heat_pump: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_min_1[s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_min[s,g,t],
                        self.variables.mu_min[s,g,t]])

                    self.constraints.SOS_min_2[s,g,t] = m.addConstr(
                        self.variables.SOS_min[s,g,t],
                        gb.GRB.EQUAL,
                        self.variables.HP_load[s,g,t],
                        name='min elec production SOS 2 ({0},{1},{2})'.format(s,g,t))     
                        
        for t in time:
            for g in CHP_sorted['ex']: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_min_1['up',s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_min[s,g,t],
                        self.variables.mu_min[s,g,t]])

                    self.constraints.SOS_min_2['up',s,g,t] = m.addConstr(
                        self.variables.SOS_min[s,g,t],
                        gb.GRB.EQUAL,
                        self.variables.P[s,g,t],
                        name='min elec production SOS 2 ({0},{1},{2})'.format(s,g,t)) 
                        

        for t in time:
            for g in CHP: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_min_1['0',s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_min_0[s,g,t],
                        self.variables.mu_min_0[s,g,t]])

                    self.constraints.SOS_min_2['0',s,g,t] = m.addConstr(
                        self.variables.SOS_min_0[s,g,t],
                        gb.GRB.EQUAL,
                        self.variables.P_0[s,g,t],
                        name='min elec production SOS 2 0 ({0},{1},{2})'.format(s,g,t)) 
                        
                        
                        
        self.constraints.SOS_max_1 = {}
        self.constraints.SOS_max_2 = {}

        for t in time:
            for g in gen: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_max_1[s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_max[s,g,t],
                        self.variables.mu_max[s,g,t]])

                    self.constraints.SOS_max_2[s,g,t] = m.addConstr(
                        self.variables.SOS_max[s,g,t],
                        gb.GRB.EQUAL,
                        -self.variables.P[s,g,t]+self.data.elec_maxprod[g],
                        name='max prod SOS 2 ({0},{1},{2})'.format(s,g,t)) 
        for t in time:
            for g in wind: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_max_1[s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_max[s,g,t],
                        self.variables.mu_max[s,g,t]])
                        
                    self.constraints.SOS_max_2[s,g,t] = m.addConstr(
                        self.variables.SOS_max[s,g,t],
                        gb.GRB.EQUAL,
                        -self.variables.P[s,g,t]+self.data.W_max*self.data.wind_scenario[s,g,t],
                        name='max prod SOS 2 ({0},{1},{2})'.format(s,g,t)) 


        for t in time:
            for g in heat_pump: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_max_1[s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_max[s,g,t],
                        self.variables.mu_max[s,g,t]])

                    self.constraints.SOS_max_2[s,g,t] = m.addConstr(
                        self.variables.SOS_max[s,g,t],
                        gb.GRB.EQUAL,
                        -self.variables.HP_load[s,g,t]+self.variables.HP_max_load[g,t],
                        name='max prod SOS 2 ({0},{1},{2})'.format(s,g,t)) 


                        
        for t in time:
            for g in CHP_sorted['ex']: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_max_1['up',s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_max[s,g,t],
                        self.variables.mu_max[s,g,t]])

                    self.constraints.SOS_max_2['up',s,g,t] = m.addConstr(
                        self.variables.SOS_max[s,g,t],
                        gb.GRB.EQUAL,
                        -self.variables.P[s,g,t]+((self.data.CHP_maxprod[g]-self.data.rho_heat[g]*self.variables.Q[g,t])/self.data.rho_elec[g]-self.data.r_min[g]*self.variables.Q[g,t]),
                        name='max prod SOS 2 ({0},{1},{2})'.format(s,g,t)) 
                        
        for t in time:
            for g in CHP: 
                for s in scenario_kmeans_0:
                    
                    self.constraints.SOS_max_1['0',s,g,t] = m.addSOS(
                        gb.GRB.SOS_TYPE1,
                        [self.variables.SOS_max_0[s,g,t],
                        self.variables.mu_max_0[s,g,t]])
                    
                    self.constraints.SOS_max_2['0',s,g,t] = m.addConstr(
                        self.variables.SOS_max_0[s,g,t],
                        gb.GRB.EQUAL,
                        -self.variables.P_0[s,g,t]+self.variables.Q[g,t]*self.data.r_min[g],
                        name='max prod SOS 2 0 ({0},{1},{2})'.format(s,g,t)) 
                        
#%% SOLVE THE MPEC (non decomposed)

Q_hierarchical_0 = {}
storage_energy_hierarchical_0 = {}
storage_charge_hierarchical_0 = {}
storage_discharge_hierarchical_0 = {}
elec_maxprod_hierarchical_0 = {}
elec_minprod_hierarchical_0 = {}

for W_max in W_range:                 
    heat_dispatch = hierarchical_heat(W_max)
    heat_dispatch.model.params.OutputFlag = 0
    heat_dispatch.optimize()

    
    for t in time:
        
        for g in CHP+heat_only:
            Q_hierarchical_0[W_max,g,t]=heat_dispatch.variables.Q[g,t].x
            
        for g in heat_pump:
            Q_hierarchical_0[W_max,g,t]=heat_dispatch.variables.HP_max_load[g,t].x*COP[g]

        for g in heat_storage:
            Q_hierarchical_0[W_max,g,t]=heat_dispatch.variables.storage_discharge[g,t].x - heat_dispatch.variables.storage_charge[g,t].x   
            storage_energy_hierarchical_0[W_max,g,t] = heat_dispatch.variables.storage_energy[g,t].x
            storage_charge_hierarchical_0[W_max,g,t] = heat_dispatch.variables.storage_charge[g,t].x
            storage_discharge_hierarchical_0[W_max,g,t] = heat_dispatch.variables.storage_discharge[g,t].x

    for t in time:
        
        for h in CHP_sorted['bp']:
            elec_maxprod_hierarchical_0[W_max,h,t] = 0
    
        for h in CHP_sorted['ex']:
            elec_maxprod_hierarchical_0[W_max,h,t] = (CHP_maxprod[h]-rho_heat[h]*Q_hierarchical_0[W_max,h,t])/rho_elec[h] - r_min[h]*Q_hierarchical_0[W_max,h,t]    

        for g in gen:
            elec_maxprod_hierarchical_0[W_max,g,t] = elec_maxprod[g]   

        for g in wind:
            elec_maxprod_hierarchical_0[W_max,g,t] = W_max 
            
        for g in heat_pump:
            elec_maxprod_hierarchical_0[W_max,g,t] = Q_hierarchical_0[W_max,g,t]/COP[g]
            
        for h in CHP:
            elec_minprod_hierarchical_0[W_max,h,t] = r_min[h]*Q_hierarchical_0[W_max,h,t]                               
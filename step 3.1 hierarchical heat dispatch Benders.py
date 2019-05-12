"""
Created on Wed Jun 22 08:48:21 2016

@author: lemitri
"""

import os
import pandas as pd
import scipy.stats as sp
#import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
sb.set_style('ticks')

#import time as ttt

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
    
    
class Benders_subproblem1:
    
    def __init__(self, MP, s0):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self.results = expando()
        
        self.MP = MP
        self.data.s0 = s0

        self._build_model()
        #self.update_constr()
        #self.update_fixed_vars()
    
    def optimize(self):
        self.model.optimize()
    
   
    def _build_model(self):
        
        self.model = gb.Model()
        #self.model.setParam('OutputFlag', False)
        self._build_variables()
        self._build_objective()
        self._build_constraints()
        self.model.update()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        heat_storage= self.MP.data.heat_storage
        heat_pump= self.MP.data.heat_pump
        m = self.model




        ## COMPLICTING VARIABLES 

        self.variables.Q = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q({0},{1})'.format(h,t))
                
        self.variables.storage_discharge = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_discharge[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage discharge({0},{1})'.format(h,t))
                    
        self.variables.storage_charge = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_charge[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage charge({0},{1})'.format(h,t))
                    

        self.variables.storage_energy = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                self.variables.storage_energy[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxcapacity[h],name='storage energy({0},{1})'.format(h,t))                


        self.variables.HP_max_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_max_load[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h]/self.MP.data.COP[h],name='heat pump load max({0},{1})'.format(h,t))                              
       
        ## 2nd stage variables

        self.variables.Q_down = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_down[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q down({0},{1})'.format(h,t))

        self.variables.Q_up = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_up[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q up({0},{1})'.format(h,t))


        self.variables.storage_discharge_scenario = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_discharge_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage discharge scenario({0},{1})'.format(h,t))
                    
        self.variables.storage_charge_scenario = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_charge_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage charge scenario({0},{1})'.format(h,t))
                    

        self.variables.storage_energy_scenario = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                self.variables.storage_energy_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxcapacity[h],name='storage energy scenario({0},{1})'.format(h,t))                

        # electricity market optimization variables

        ## primal variables

        self.variables.HP_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_load[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h]/self.MP.data.COP[h],name='heat pump load max({0},{1})'.format(h,t))                


        self.variables.P = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in gen+wind+CHP_sorted['ex']:
                self.variables.P[g,t] = m.addVar(lb=0) # dispatch of electricity generators

        self.variables.P_0 = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP:
                self.variables.P_0[g,t] = m.addVar(lb=0) # dispatch of electricity generators
                  
                   
        ## dual variables
        
        self.variables.lambda_powerbalance = {}
        for t in time:
            self.variables.lambda_powerbalance[t] = m.addVar(lb=-gb.GRB.INFINITY,name='lambda powerbalance({0})'.format(t))
        
        self.variables.mu_max = {}
        for t in time:
            for g in CHP_sorted['ex']+gen+wind+heat_pump:
                self.variables.mu_max[g,t] = m.addVar(lb=0,name='mu max({0},{1})'.format(g,t)) 


        self.variables.mu_max_0 = {}
        for t in time:
            for g in CHP:
                self.variables.mu_max_0[g,t] = m.addVar(lb=0,name='mu max 0({0},{1})'.format(g,t)) 
        
        
        m.update()
        
    def _build_objective(self): # building the objective function for the heat maret clearing
        
        #indexes shortcuts 
        #S_all = self.MP.data.S_all
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        s0 = self.data.s0
        m = self.model
        
        #Objective function


        m.setObjective(gb.quicksum(self.MP.data.alpha_up[s0,h,t]*self.variables.Q_up[h,t] - self.MP.data.alpha_down[s0,h,t]*self.variables.Q_down[h,t] for h in CHP+heat_only for t in time) + gb.quicksum(self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]*self.variables.P[h,t] for h in CHP for t in time) + gb.quicksum(self.MP.data.alpha[s0,h,t]*self.variables.P[h,t] for h in gen+wind for t in time) + gb.quicksum(self.variables.mu_max[h,t]*self.MP.data.elec_maxprod[h] for h in gen for t in time) + gb.quicksum(self.variables.mu_max[h,t]*self.MP.data.wind_scenario[s0,h,t]*self.MP.data.W_max for h in wind for t in time) - gb.quicksum(self.variables.lambda_powerbalance[t]*self.MP.data.elec_load[s0,t] for t in time),   
            gb.GRB.MINIMIZE)            
        
    def _build_constraints(self):
        
        #indexes shortcuts 
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        heat_storage= self.MP.data.heat_storage
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

       #2nd stage constraints

        #CHP's joint FOR in each scenario 
        
        self.constraints.heat_minprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_minprod[h,t] = m.addConstr(
                    self.variables.Q[h,t] + (self.variables.Q_up[h,t]-self.variables.Q_down[h,t]),
                    gb.GRB.GREATER_EQUAL,
                    0) 

        self.constraints.heat_maxprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_maxprod[h,t] = m.addConstr(
                    self.variables.Q[h,t]+(self.variables.Q_up[h,t]-self.variables.Q_down[h,t]),
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.heat_maxprod[h])
                        
        self.constraints.CHP_maxprod = {} 
        self.constraints.CHP_ratio = {}
        
        for t in time:
            for h in CHP_sorted['ex']:
                
                self.constraints.CHP_maxprod[h,t] = m.addConstr(
                    (self.variables.Q[h,t] + self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.rho_heat[h]+self.MP.data.rho_elec[h]*(self.variables.P[h,t]+self.variables.P_0[h,t]),
                     gb.GRB.LESS_EQUAL,
                     self.MP.data.CHP_maxprod[h])
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    (self.variables.P[h,t]+self.variables.P_0[h,t]),
                    gb.GRB.GREATER_EQUAL,
                    (self.variables.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.r_min[h])

            for h in CHP_sorted['bp']:
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    self.variables.P_0[h,t],
                    gb.GRB.EQUAL,
                    (self.variables.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.r_min[h])


        self.constraints.heat_balance_scenario = {} 
        
        for t in time:
            self.constraints.heat_balance_scenario[t] = m.addConstr(
                    gb.quicksum(self.variables.Q[h,t] + self.variables.Q_up[h,t] - self.variables.Q_down[h,t] for h in CHP+heat_only) + gb.quicksum(self.MP.data.COP[h]*self.variables.HP_load[h,t] for h in heat_pump) + gb.quicksum(self.variables.storage_discharge_scenario[h,t]-self.variables.storage_charge_scenario[h,t] for h in heat_storage),
                    gb.GRB.EQUAL,
                    self.MP.data.heat_load[t])
                    

        # storage (and stage)

        # heat storage: storage states update

        self.constraints.storage_update_scenario={}
        self.constraints.storage_init_scenario={}
        self.constraints.storage_final_scenario={}
        
        for h in heat_storage:        
            for (t1,t2) in zip(time[:-1],time[1:]):
                self.constraints.storage_update_scenario[h,t2]=m.addConstr(
                    self.variables.storage_energy_scenario[h,t2],
                    gb.GRB.EQUAL,
                    self.variables.storage_energy_scenario[h,t1]-self.MP.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[h,t2]+self.MP.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[h,t2]-self.MP.data.storage_loss[h])

            self.constraints.storage_init_scenario[h]=m.addConstr(
                self.variables.storage_energy_scenario[h,time[0]],
                gb.GRB.EQUAL,
                self.MP.data.storage_energy_init[h]-self.MP.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[h,time[0]]+self.MP.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[h,time[0]]-self.MP.data.storage_loss[h])
   
            self.constraints.storage_final_scenario[h]=m.addConstr(
                self.variables.storage_energy_scenario[h,time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.MP.data.storage_energy_init[h])
                    
        # Lower level Constraints
        
        ## Equality
        
        self.constraints.elec_balance = {} 
        
        for t in time:
            self.constraints.elec_balance[t] = m.addConstr(
                gb.quicksum(self.variables.P[g,t] for g in CHP_sorted['ex']+gen+wind)+gb.quicksum(self.variables.P_0[g,t] for g in CHP)-gb.quicksum(self.variables.HP_load[h,t] for h in heat_pump),
                gb.GRB.EQUAL,
                self.MP.data.elec_load[s0,t])     
                
        ## Inequality
        
        self.constraints.P_max = {}
        
        for t in time:  


            for g in gen:
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.elec_maxprod[g])
       
            for g in wind:                    
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.wind_scenario[s0,g,t]*self.MP.data.W_max)

            for g in CHP_sorted['ex']:
                self.constraints.P_max['up',g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    (self.MP.data.CHP_maxprod[g]-self.MP.data.rho_heat[g]*self.variables.Q[g,t])/self.MP.data.rho_elec[g]-self.variables.Q[g,t]*self.MP.data.r_min[g])                      
                                            
                    
            for g in CHP:
                self.constraints.P_max['0',g,t] = m.addConstr(
                    self.variables.P_0[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.variables.Q[g,t]*self.MP.data.r_min[g]) 


            for g in heat_pump:
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.HP_load[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.variables.HP_max_load[g,t])
                        
        # lower_level_problem: dual feasibility of lower level problem: Delta L =0 (equivalent to dual constraints)
             
        self.constraints.dL = {} # gradient of lagragian function / dP_gen[g,t]
        
        for t in time:
                
            for g in gen+wind:
                self.constraints.dL[g,t] = m.addConstr(
                    self.MP.data.alpha[s0,g,t]-self.variables.lambda_powerbalance[t]+self.variables.mu_max[g,t],
                    gb.GRB.GREATER_EQUAL,
                    0)
                                        
            for h in CHP:
                self.constraints.dL['0',h,t] = m.addConstr(
                    self.MP.data.alpha_min-self.variables.lambda_powerbalance[t]+self.variables.mu_max_0[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

            for h in CHP_sorted['ex']:
                self.constraints.dL['up',h,t] = m.addConstr(
                    self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]-self.variables.lambda_powerbalance[t]+self.variables.mu_max[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)
                                         

            for h in heat_pump:
                self.constraints.dL[h,t] = m.addConstr(
                    -self.MP.data.alpha_HP+self.variables.lambda_powerbalance[t]+self.variables.mu_max[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

        
        a = -1
        b = -1   
        c = -1  
        d = 1
                                                      
        self.constraints.strong_duality=m.addConstr(
            gb.quicksum(-self.variables.mu_max[g,t]*self.MP.data.elec_maxprod[g] for g in gen for t in time) + gb.quicksum(-self.variables.mu_max[g,t]*self.MP.data.wind_scenario[s0,g,t]*self.MP.data.W_max for g in wind for t in time)
            + gb.quicksum(self.variables.mu_max[h,t]*a for h in CHP_sorted['ex'] for t in time) + gb.quicksum(self.variables.mu_max_0[h,t]*b for h in CHP for t in time) + gb.quicksum(self.variables.mu_max[h,t]*c for h in heat_pump for t in time) + gb.quicksum(self.variables.lambda_powerbalance[t]*self.MP.data.elec_load[s0,t] for t in time)
            + gb.quicksum(-self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]*self.variables.P[h,t] for h in CHP_sorted['ex'] for t in time) + gb.quicksum(-self.MP.data.alpha_min*self.variables.P_0[h,t] for h in CHP for t in time) + gb.quicksum(-self.MP.data.alpha[s0,g,t]*self.variables.P[g,t] for g in gen+wind for t in time) + gb.quicksum(self.variables.HP_load[h,t]*self.MP.data.alpha_HP for h in heat_pump for t in time),
            gb.GRB.EQUAL,
            0,
            name='strong_duality1')
            
            # parameters a b c and d will be fixed to:
            # d = + self.MP.variables.alpha_HP[h,t].x
            # c = - self.MP.variables.HP_max_load[h,t].x 
            # b = - self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x
            # a = - ((self.MP.data.CHP_maxprod[h]-self.MP.variables.Q[h,t].x*self.MP.data.rho_heat[h])/self.MP.data.rho_elec[h]-self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x)


        self.constraints.fix_Q = {}

        for t in time:
            for h in CHP+heat_only:
                self.constraints.fix_Q[h,t] = m.addConstr(self.variables.Q[h,t],gb.GRB.EQUAL,0)

        self.constraints.fix_HP_max_load = {}

        for t in time:        
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t] = m.addConstr(self.variables.HP_max_load[h,t],gb.GRB.EQUAL,0)

        self.constraints.fix_storage_energy = {}

        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t] = m.addConstr(self.variables.storage_energy[h,t],gb.GRB.EQUAL,0) 

        self.constraints.fix_storage_discharge = {}
        
        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_discharge[h,t] = m.addConstr(self.variables.storage_discharge[h,t],gb.GRB.EQUAL,0) 

        self.constraints.fix_storage_charge = {}
        
        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_charge[h,t] = m.addConstr(self.variables.storage_charge[h,t],gb.GRB.EQUAL,0) 
                
                
        m.update()
            

    def update_fixed_vars_warm_start(self): 

        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        heat_only=self.MP.data.heat_only
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

        
        for t in time:
            for h in CHP_sorted['ex']:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max[h,t],-((self.MP.data.CHP_maxprod[h]-self.MP.data.Q_hierarchical_0[h,t]*self.MP.data.rho_heat[h])/self.MP.data.rho_elec[h]-self.MP.data.r_min[h]*self.MP.data.Q_hierarchical_0[h,t])) 
            for h in CHP:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max_0[h,t],-self.MP.data.r_min[h]*self.MP.data.Q_hierarchical_0[h,t]) 
            for h in heat_pump:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max[h,t],-self.MP.data.Q_hierarchical_0[h,t]/self.MP.data.COP[h])
                
            # d = + self.MP.variables.alpha_HP[h,t].x
            # c = - self.MP.variables.HP_max_load[h,t].x 
            # b = - self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x
            # a = - ((self.MP.data.CHP_maxprod[h]-self.MP.variables.Q[h,t].x*self.MP.data.rho_heat[h])/self.MP.data.rho_elec[h]-self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x)

            
                
        for t in time:
                
            for h in CHP:                
                self.constraints.fix_Q[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]

            for h in heat_only:
                self.constraints.fix_Q[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]

      
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]/self.MP.data.COP[h]
                
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t].rhs = self.MP.data.storage_energy_hierarchical_0[h,t]
                self.constraints.fix_storage_discharge[h,t].rhs = self.MP.data.storage_discharge_hierarchical_0[h,t]
                self.constraints.fix_storage_charge[h,t].rhs = self.MP.data.storage_charge_hierarchical_0[h,t]

        m.update()
        
        
    def update_fixed_vars(self): 

        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        heat_only=self.MP.data.heat_only
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

        
        for t in time:
            for h in CHP_sorted['ex']:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max[h,t],-((self.MP.data.CHP_maxprod[h]-self.MP.variables.Q[h,t].x*self.MP.data.rho_heat[h])/self.MP.data.rho_elec[h]-self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x)) 
            for h in CHP:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max_0[h,t],-self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x) 
            for h in heat_pump:
                m.chgCoeff(self.constraints.strong_duality,self.variables.mu_max[h,t],-self.MP.variables.HP_max_load[h,t].x)
                
            # d = + self.MP.variables.alpha_HP[h,t].x
            # c = - self.MP.variables.HP_max_load[h,t].x 
            # b = - self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x
            # a = - ((self.MP.data.CHP_maxprod[h]-self.MP.variables.Q[h,t].x*self.MP.data.rho_heat[h])/self.MP.data.rho_elec[h]-self.MP.data.r_min[h]*self.MP.variables.Q[h,t].x)

            
                
        for t in time:
                
            for h in CHP:                
                self.constraints.fix_Q[h,t].rhs = self.MP.variables.Q[h,t].x

            for h in heat_only:
                self.constraints.fix_Q[h,t].rhs = self.MP.variables.Q[h,t].x

      
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t].rhs = self.MP.variables.HP_max_load[h,t].x
                
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t].rhs = self.MP.variables.storage_energy[h,t].x
                self.constraints.fix_storage_discharge[h,t].rhs = self.MP.variables.storage_discharge[h,t].x
                self.constraints.fix_storage_charge[h,t].rhs = self.MP.variables.storage_charge[h,t].x

        m.update()

#%%

class Benders_subproblem2:
    
    def __init__(self, MP, s0):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self.results = expando()

        self.data.s0 = s0        
        self.MP = MP

        self.SUB1 = self.MP.submodels1[s0]
        
        self._build_model()


    def optimize(self):
        self.model.optimize()
    
   
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
        self.model.update()
    
    def _build_variables(self):
        
        #indexes shortcuts 
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        heat_storage= self.MP.data.heat_storage
        heat_pump= self.MP.data.heat_pump
        m = self.model

        ## COMPLICTING VARIABLES 

        self.variables.Q = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q({0},{1})'.format(h,t))
                
        self.variables.storage_discharge = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_discharge[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage discharge({0},{1})'.format(h,t))
                    
        self.variables.storage_charge = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_charge[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage charge({0},{1})'.format(h,t))
                    

        self.variables.storage_energy = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                self.variables.storage_energy[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxcapacity[h],name='storage energy({0},{1})'.format(h,t))                


        self.variables.HP_max_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_max_load[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h]/self.MP.data.COP[h],name='heat pump load max({0},{1})'.format(h,t))                
             
       
        ## 2nd stage variables

        self.variables.Q_down = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_down[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q down({0},{1})'.format(h,t))

        self.variables.Q_up = {} #heat production from CHPs and HO units (first satge)
        for t in time:
            for h in CHP+heat_only:
                self.variables.Q_up[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h],name='Q up({0},{1})'.format(h,t))


        self.variables.storage_discharge_scenario = {} #heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_discharge_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage discharge scenario({0},{1})'.format(h,t))
                    
        self.variables.storage_charge_scenario = {} #heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage:
                self.variables.storage_charge_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxprod[h],name='storage charge scenario({0},{1})'.format(h,t))
                    

        self.variables.storage_energy_scenario = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage:
                self.variables.storage_energy_scenario[h,t] = m.addVar(lb=0,ub=self.MP.data.storage_maxcapacity[h],name='storage energy scenario({0},{1})'.format(h,t))                

        # electricity market optimization variables

        ## primal variables

        self.variables.HP_load = {} #heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_pump:
                self.variables.HP_load[h,t] = m.addVar(lb=0,ub=self.MP.data.heat_maxprod[h]/self.MP.data.COP[h],name='heat pump load max({0},{1})'.format(h,t))                


        self.variables.P = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in gen+wind+CHP_sorted['ex']:
                self.variables.P[g,t] = m.addVar(lb=0) # dispatch of electricity generators

        self.variables.P_0 = {} # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP:
                self.variables.P_0[g,t] = m.addVar(lb=0) # dispatch of electricity generators
                  
                   
        ## dual variables
        
        self.variables.lambda_powerbalance = {}
        for t in time:
            self.variables.lambda_powerbalance[t] = m.addVar(lb=-gb.GRB.INFINITY,name='lambda powerbalance({0})'.format(t))
        
        self.variables.mu_max = {}
        for t in time:
            for g in CHP_sorted['ex']+gen+wind+heat_pump:
                self.variables.mu_max[g,t] = m.addVar(lb=0,name='mu max({0},{1})'.format(g,t)) 


        self.variables.mu_max_0 = {}
        for t in time:
            for g in CHP:
                self.variables.mu_max_0[g,t] = m.addVar(lb=0,name='mu max 0({0},{1})'.format(g,t)) 
        
        self.variables.r1 = m.addVar()
        self.variables.r2 = m.addVar() 
               
        m.update()                            
                
    def _build_objective(self): # building the objective function for the heat maret clearing
        
        #indexes shortcuts 
        #S_all = self.MP.data.S_all
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        s0 = self.data.s0
        m = self.model

        m.setObjective(100*(self.variables.r2+self.variables.r1)+gb.quicksum(self.MP.data.alpha_up[s0,h,t]*self.variables.Q_up[h,t] - self.MP.data.alpha_down[s0,h,t]*self.variables.Q_down[h,t] for h in CHP+heat_only for t in time) + gb.quicksum(self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]*self.variables.P[h,t] for h in CHP for t in time) + gb.quicksum(self.MP.data.alpha[s0,h,t]*self.variables.P[h,t] for h in gen+wind for t in time) + gb.quicksum(self.variables.mu_max[h,t]*self.MP.data.elec_maxprod[h] for h in gen for t in time) + gb.quicksum(self.variables.mu_max[h,t]*self.MP.data.wind_scenario[s0,h,t]*self.MP.data.W_max for h in wind for t in time) - gb.quicksum(self.variables.lambda_powerbalance[t]*self.MP.data.elec_load[s0,t] for t in time),   
            gb.GRB.MINIMIZE) 
        
    def _build_constraints(self):
        
        #indexes shortcuts 
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        heat_storage= self.MP.data.heat_storage
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

    def _build_constraints(self):
        
        #indexes shortcuts 
        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        gen=self.MP.data.gen
        wind=self.MP.data.wind
        heat_only=self.MP.data.heat_only
        heat_storage= self.MP.data.heat_storage
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

       #2nd stage constraints

        #CHP's joint FOR in each scenario 
        
        self.constraints.heat_minprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_minprod[h,t] = m.addConstr(
                    self.variables.Q[h,t] + (self.variables.Q_up[h,t]-self.variables.Q_down[h,t]),
                    gb.GRB.GREATER_EQUAL,
                    0) 

        self.constraints.heat_maxprod = {} 

        for t in time:
            for h in CHP+heat_only:
                
                self.constraints.heat_maxprod[h,t] = m.addConstr(
                    self.variables.Q[h,t]+(self.variables.Q_up[h,t]-self.variables.Q_down[h,t]),
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.heat_maxprod[h])
                        
        self.constraints.CHP_maxprod = {} 
        self.constraints.CHP_ratio = {}
        
        for t in time:
            for h in CHP_sorted['ex']:
                
                self.constraints.CHP_maxprod[h,t] = m.addConstr(
                    (self.variables.Q[h,t] + self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.rho_heat[h]+self.MP.data.rho_elec[h]*(self.variables.P[h,t]+self.variables.P_0[h,t]),
                     gb.GRB.LESS_EQUAL,
                     self.MP.data.CHP_maxprod[h])
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    (self.variables.P[h,t]+self.variables.P_0[h,t]),
                    gb.GRB.GREATER_EQUAL,
                    (self.variables.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.r_min[h])

            for h in CHP_sorted['bp']:
                
                self.constraints.CHP_ratio[h,t] = m.addConstr(
                    self.variables.P_0[h,t],
                    gb.GRB.EQUAL,
                    (self.variables.Q[h,t]+self.variables.Q_up[h,t]-self.variables.Q_down[h,t])*self.MP.data.r_min[h])


        self.constraints.heat_balance_scenario = {} 
        
        for t in time:
            self.constraints.heat_balance_scenario[t] = m.addConstr(
                    gb.quicksum(self.variables.Q[h,t] + self.variables.Q_up[h,t] - self.variables.Q_down[h,t] for h in CHP+heat_only) + gb.quicksum(self.MP.data.COP[h]*self.variables.HP_load[h,t] for h in heat_pump) + gb.quicksum(self.variables.storage_discharge_scenario[h,t]-self.variables.storage_charge_scenario[h,t] for h in heat_storage),
                    gb.GRB.EQUAL,
                    self.MP.data.heat_load[t])
                    

        # storage (and stage)

        # heat storage: storage states update

        self.constraints.storage_update_scenario={}
        self.constraints.storage_init_scenario={}
        self.constraints.storage_final_scenario={}
        
        for h in heat_storage:        
            for (t1,t2) in zip(time[:-1],time[1:]):
                self.constraints.storage_update_scenario[h,t2]=m.addConstr(
                    self.variables.storage_energy_scenario[h,t2],
                    gb.GRB.EQUAL,
                    self.variables.storage_energy_scenario[h,t1]-self.MP.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[h,t2]+self.MP.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[h,t2]-self.MP.data.storage_loss[h])

            self.constraints.storage_init_scenario[h]=m.addConstr(
                self.variables.storage_energy_scenario[h,time[0]],
                gb.GRB.EQUAL,
                self.MP.data.storage_energy_init[h]-self.MP.data.storage_discharge_eff[h]*self.variables.storage_discharge_scenario[h,time[0]]+self.MP.data.storage_charge_eff[h]*self.variables.storage_charge_scenario[h,time[0]]-self.MP.data.storage_loss[h])
   
            self.constraints.storage_final_scenario[h]=m.addConstr(
                self.variables.storage_energy_scenario[h,time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.MP.data.storage_energy_init[h])
                    
        # Lower level Constraints
        
        ## Equality
        
        self.constraints.elec_balance = {} 
        
        for t in time:
            self.constraints.elec_balance[t] = m.addConstr(
                gb.quicksum(self.variables.P[g,t] for g in CHP_sorted['ex']+gen+wind)+gb.quicksum(self.variables.P_0[g,t] for g in CHP)-gb.quicksum(self.variables.HP_load[h,t] for h in heat_pump),
                gb.GRB.EQUAL,
                self.MP.data.elec_load[s0,t])     
                
        ## Inequality
        
        self.constraints.P_max = {}
        
        for t in time:  


            for g in gen:
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.elec_maxprod[g])
       
            for g in wind:                    
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.MP.data.wind_scenario[s0,g,t]*self.MP.data.W_max)

            for g in CHP_sorted['ex']:
                self.constraints.P_max['up',g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    (self.MP.data.CHP_maxprod[g]-self.MP.data.rho_heat[g]*self.variables.Q[g,t])/self.MP.data.rho_elec[g]-self.variables.Q[g,t]*self.MP.data.r_min[g])                      
                                            
                    
            for g in CHP:
                self.constraints.P_max['0',g,t] = m.addConstr(
                    self.variables.P_0[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.variables.Q[g,t]*self.MP.data.r_min[g]) 


            for g in heat_pump:
                self.constraints.P_max[g,t] = m.addConstr(
                    self.variables.HP_load[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.variables.HP_max_load[g,t])
                        
        # lower_level_problem: dual feasibility of lower level problem: Delta L =0 (equivalent to dual constraints)
             
        self.constraints.dL = {} # gradient of lagragian function / dP_gen[g,t]
        
        for t in time:
                
            for g in gen+wind:
                self.constraints.dL[g,t] = m.addConstr(
                    self.MP.data.alpha[s0,g,t]-self.variables.lambda_powerbalance[t]+self.variables.mu_max[g,t],
                    gb.GRB.GREATER_EQUAL,
                    0)
                                        
            for h in CHP:
                self.constraints.dL['0',h,t] = m.addConstr(
                    self.MP.data.alpha_min-self.variables.lambda_powerbalance[t]+self.variables.mu_max_0[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

            for h in CHP_sorted['ex']:
                self.constraints.dL['up',h,t] = m.addConstr(
                    self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]-self.variables.lambda_powerbalance[t]+self.variables.mu_max[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)
                                         

            for h in heat_pump:
                self.constraints.dL[h,t] = m.addConstr(
                    -self.MP.data.alpha_HP+self.variables.lambda_powerbalance[t]+self.variables.mu_max[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

        
        a = 1
        b = -1   
        c = -1  
        d = 1
                                                      
        self.constraints.strong_duality=m.addConstr(
            (-self.variables.r2+self.variables.r1)+gb.quicksum(-self.variables.mu_max[g,t]*self.MP.data.elec_maxprod[g] for g in gen for t in time) + gb.quicksum(-self.variables.mu_max[g,t]*self.MP.data.wind_scenario[s0,g,t]*self.MP.data.W_max for g in wind for t in time)
            + gb.quicksum(-self.variables.mu_max[h,t]*self.MP.data.CHP_maxprod[h]/self.MP.data.rho_elec[h] for h in CHP_sorted['ex'] for t in time) + gb.quicksum(self.variables.Q[h,t]*a for h in CHP_sorted['ex'] for t in time) + gb.quicksum(self.variables.Q[h,t]*b for h in CHP for t in time) + gb.quicksum(self.variables.HP_max_load[h,t]*c for h in heat_pump for t in time) + gb.quicksum(self.variables.lambda_powerbalance[t]*self.MP.data.elec_load[s0,t] for t in time)
            + gb.quicksum(-self.MP.data.alpha[s0,h,t]*self.MP.data.rho_elec[h]*self.variables.P[h,t] for h in CHP_sorted['ex'] for t in time) + gb.quicksum(-self.MP.data.alpha_min*self.variables.P_0[h,t] for h in CHP for t in time) + gb.quicksum(-self.MP.data.alpha[s0,g,t]*self.variables.P[g,t] for g in gen+wind for t in time) + gb.quicksum(self.variables.HP_load[h,t]*self.MP.data.alpha_HP for h in heat_pump for t in time),
            gb.GRB.EQUAL,
            0,
            name='strong_duality1') # how to update only the parameters??????
            # d = + self.SUB1.variables.HP_load[h,t].x
            # c = - self.SUB1.variables.mu_max[h,t].x 
            # b = - self.MP.data.r_min[h]*self.SUB1.variables.mu_max_0[h,t].x
            # a = + self.SUB1.variables.mu_max[h,t].x*(self.MP.data.rho_heat[h]/self.MP.data.rho_elec[h]+self.MP.data.r_min[h])



        self.constraints.fix_Q = {}

        for t in time:
            for h in CHP+heat_only:
                self.constraints.fix_Q[h,t] = m.addConstr(self.variables.Q[h,t],gb.GRB.EQUAL,0)

        self.constraints.fix_HP_max_load = {}

        for t in time:        
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t] = m.addConstr(self.variables.HP_max_load[h,t],gb.GRB.EQUAL,0)

        self.constraints.fix_storage_energy = {}

        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t] = m.addConstr(self.variables.storage_energy[h,t],gb.GRB.EQUAL,0) 

        self.constraints.fix_storage_discharge = {}
        
        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_discharge[h,t] = m.addConstr(self.variables.storage_discharge[h,t],gb.GRB.EQUAL,0) 

        self.constraints.fix_storage_charge = {}
        
        for t in time:        
            for h in heat_storage:
                self.constraints.fix_storage_charge[h,t] = m.addConstr(self.variables.storage_charge[h,t],gb.GRB.EQUAL,0) 

        # BECAUSE WE LINEARIZE THE STRONG DUALITY

        self.constraints.fix_HP_load = {}

        for t in time:
            for h in heat_pump:
                self.constraints.fix_HP_load[h,t] = m.addConstr(self.variables.HP_load[h,t],gb.GRB.EQUAL,0)
                

        self.constraints.fix_mu_max = {}

        for t in time:
            for h in CHP_sorted['ex']+heat_pump:
                self.constraints.fix_mu_max[h,t] = m.addConstr(self.variables.mu_max[h,t],gb.GRB.EQUAL,0)

        self.constraints.fix_mu_max_0 = {}

        for t in time:
            for h in CHP:
                self.constraints.fix_mu_max_0[h,t] = m.addConstr(self.variables.mu_max_0[h,t],gb.GRB.EQUAL,0)
                
        m.update()
            


    def update_fixed_vars_warm_start(self): 

        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        heat_only=self.MP.data.heat_only
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

        
        for t in time:
            for h in CHP_sorted['ex']:
                m.chgCoeff(self.constraints.strong_duality,self.variables.Q[h,t],self.SUB1.variables.mu_max[h,t].x*(self.MP.data.rho_heat[h]/self.MP.data.rho_elec[h]+self.MP.data.r_min[h])) 
            for h in CHP:
                m.chgCoeff(self.constraints.strong_duality,self.variables.Q[h,t],- self.MP.data.r_min[h]*self.SUB1.variables.mu_max_0[h,t].x) 
            for h in heat_pump:
                m.chgCoeff(self.constraints.strong_duality,self.variables.HP_max_load[h,t],- self.SUB1.variables.mu_max[h,t].x)

            
                
        for t in time:
                
            for h in CHP:                
                self.constraints.fix_Q[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]

            for h in heat_only:
                self.constraints.fix_Q[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]

      
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t].rhs = self.MP.data.Q_hierarchical_0[h,t]/self.MP.data.COP[h]
                
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t].rhs = self.MP.data.storage_energy_hierarchical_0[h,t]
                self.constraints.fix_storage_discharge[h,t].rhs = self.MP.data.storage_discharge_hierarchical_0[h,t]
                self.constraints.fix_storage_charge[h,t].rhs = self.MP.data.storage_charge_hierarchical_0[h,t]



        for t in time:
            
            for h in heat_pump:
                self.constraints.fix_HP_load[h,t].rhs = self.SUB1.variables.HP_load[h,t].x

            for h in CHP_sorted['ex']+heat_pump:
                self.constraints.fix_mu_max[h,t].rhs = self.SUB1.variables.mu_max[h,t].x

            for h in CHP:
                self.constraints.fix_mu_max_0[h,t].rhs = self.SUB1.variables.mu_max_0[h,t].x  
                
                
                
        m.update()
        
        
        
        
    def update_fixed_vars(self): 

        time = self.MP.data.time
        CHP = self.MP.data.CHP
        CHP_sorted = self.MP.data.CHP_sorted
        heat_only=self.MP.data.heat_only
        heat_pump= self.MP.data.heat_pump
        s0 = self.data.s0
        m = self.model

        
        for t in time:
            for h in CHP_sorted['ex']:
                m.chgCoeff(self.constraints.strong_duality,self.variables.Q[h,t],self.SUB1.variables.mu_max[h,t].x*(self.MP.data.rho_heat[h]/self.MP.data.rho_elec[h]+self.MP.data.r_min[h])) 
            for h in CHP:
                m.chgCoeff(self.constraints.strong_duality,self.variables.Q[h,t],- self.MP.data.r_min[h]*self.SUB1.variables.mu_max_0[h,t].x) 
            for h in heat_pump:
                m.chgCoeff(self.constraints.strong_duality,self.variables.HP_max_load[h,t],- self.SUB1.variables.mu_max[h,t].x)

            # d = + self.SUB1.variables.HP_load[h,t].x
            # c = - self.SUB1.variables.mu_max[h,t].x 
            # b = - self.MP.data.r_min[h]*self.SUB1.variables.mu_max_0[h,t].x
            # a = + self.SUB1.variables.mu_max[h,t].x*(self.MP.data.rho_heat[h]/self.MP.data.rho_elec[h]+self.MP.data.r_min[h])



                
        for t in time:
                
            for h in CHP:                
                self.constraints.fix_Q[h,t].rhs = self.MP.variables.Q[h,t].x

            for h in heat_only:
                self.constraints.fix_Q[h,t].rhs = self.MP.variables.Q[h,t].x

      
            for h in heat_pump:
                self.constraints.fix_HP_max_load[h,t].rhs = self.MP.variables.HP_max_load[h,t].x
                
            for h in heat_storage:
                self.constraints.fix_storage_energy[h,t].rhs = self.MP.variables.storage_energy[h,t].x
                self.constraints.fix_storage_discharge[h,t].rhs = self.MP.variables.storage_discharge[h,t].x
                self.constraints.fix_storage_charge[h,t].rhs = self.MP.variables.storage_charge[h,t].x


        for t in time:
            
            for h in heat_pump:
                self.constraints.fix_HP_load[h,t].rhs = self.SUB1.variables.HP_load[h,t].x

            for h in CHP_sorted['ex']+heat_pump:
                self.constraints.fix_mu_max[h,t].rhs = self.SUB1.variables.mu_max[h,t].x

            for h in CHP:
                self.constraints.fix_mu_max_0[h,t].rhs = self.SUB1.variables.mu_max_0[h,t].x  
                
        m.update()


#%%

class Benders_master:
    def __init__(self,W_max):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self.results = expando()
        self._load_data(W_max)
        self._init_benders_params(epsilon=1,max_iters=300)
        self._build_model()
    
    def master_optimize(self):
        
        time = self.data.time
        CHP = self.data.CHP
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        m = self.model

        scenario_kmeans=self.data.scenario_kmeans
        S_all_kmeans = self.data.S_all_kmeans
        
        # initial solutions: warm start

        self.data.Q_ref = self.data.Q_hierarchical_0
        self.data.HP_max_load_ref = self.data.HP_max_load_hierarchical_0
        self.data.storage_discharge_ref = self.data.storage_discharge_hierarchical_0
        self.data.storage_charge_ref = self.data.storage_charge_hierarchical_0
        self.data.storage_energy_ref = self.data.storage_energy_hierarchical_0
        
        
        # subproblems1: BUILD
        
        self.submodels1 = {s:Benders_subproblem1(self,s0=s) for s in self.data.scenario_kmeans}

        for s in scenario_kmeans:
            self.submodels1[s].model.params.OutputFlag = 0         
       
        # Update fixed variables for submodels 1 and solve.
        [self.submodels1[s].update_fixed_vars_warm_start() for s in self.data.scenario_kmeans]
        [self.submodels1[s].optimize() for s in self.data.scenario_kmeans]
        
        #subproblems 2: BUILD
        
        self.submodels2 = {s:Benders_subproblem2(self,s0=s) for s in self.data.scenario_kmeans}
        '''
        for s in self.data.scenario_kmeans:        
            self.submodels2[s].model.params.FeasibilityTol=1e-2
            self.submodels2[s].model.params.OptimalityTol=1e-2
        '''
        for s in scenario_kmeans:        
            self.submodels2[s].model.params.OutputFlag = 0 
            
        # Update fixed variables for submodels 2 and solve.            
        [self.submodels2[s].update_fixed_vars_warm_start() for s in self.data.scenario_kmeans]
        [self.submodels2[s].optimize() for s in self.data.scenario_kmeans]
        
        self.data.level = 10
        
        self.model.setObjective(gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_heat[h]*self.variables.Q[h,t] for s in scenario_kmeans for h in CHP+heat_only for t in time)/S_all_kmeans + self.variables.beta + 1/2*self.data.level*(gb.quicksum((self.variables.Q[h,t]-self.data.Q_ref[h,t])*(self.variables.Q[h,t]-self.data.Q_ref[h,t]) for h in CHP+heat_only for t in time) + gb.quicksum((self.variables.HP_max_load[h,t]-self.data.HP_max_load_ref[h,t])*(self.variables.HP_max_load[h,t]-self.data.HP_max_load_ref[h,t]) for h in heat_pump for t in time) + gb.quicksum((self.variables.storage_energy[h,t]-self.data.storage_energy_ref[h,t])*(self.variables.storage_energy[h,t]-self.data.storage_energy_ref[h,t]) + (self.variables.storage_charge[h,t]-self.data.storage_charge_ref[h,t])*(self.variables.storage_charge[h,t]-self.data.storage_charge_ref[h,t]) + (self.variables.storage_discharge[h,t]-self.data.storage_discharge_ref[h,t])*(self.variables.storage_discharge[h,t]-self.data.storage_discharge_ref[h,t]) for h in heat_storage for t in time) ),   
            gb.GRB.MINIMIZE)
            
        cut = len(self.data.cutlist)
        self.data.cutlist.append(cut)

        # Get sensitivities from subproblems 2
        self.data.sens_Q = {(h,t): sum([self.submodels2[s].constraints.fix_Q[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in CHP+heat_only for t in time}
        self.data.sens_storage_discharge = {(h,t): sum([self.submodels2[s].constraints.fix_storage_discharge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_charge = {(h,t): sum([self.submodels2[s].constraints.fix_storage_charge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_energy = {(h,t): sum([self.submodels2[s].constraints.fix_storage_energy[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_HP_max_load  =  {(h,t): sum([self.submodels2[s].constraints.fix_HP_max_load[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_pump for t in time}


        self.data.sens_Q_ref = self.data.sens_Q
        self.data.sens_storage_discharge_ref = self.data.sens_storage_discharge
        self.data.sens_storage_charge_ref = self.data.sens_storage_charge
        self.data.sens_storage_energy_ref = self.data.sens_storage_energy
        self.data.sens_HP_max_load_ref  = self.data.sens_HP_max_load


        self.data.z_sub = sum([self.submodels1[s].model.ObjVal for s in scenario_kmeans])/S_all_kmeans
        
        self.constraints.cuts[cut] = self.model.addConstr(
            self.variables.beta,
            gb.GRB.GREATER_EQUAL,
            self.data.z_sub + gb.quicksum(self.data.sens_Q[h,t]*(self.variables.Q[h,t] - self.data.Q_hierarchical_0[h,t]) for h in CHP+heat_only for t in time) + gb.quicksum(self.data.sens_storage_discharge[h,t]*(self.variables.storage_discharge[h,t] - self.data.storage_discharge_hierarchical_0[h,t] ) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_storage_charge[h,t]*(self.variables.storage_charge[h,t] - self.data.storage_charge_hierarchical_0[h,t]) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_storage_energy[h,t]*(self.variables.storage_energy[h,t] - self.data.storage_energy_hierarchical_0[h,t] ) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_HP_max_load[h,t]*(self.variables.HP_max_load[h,t] - self.data.Q_hierarchical_0[h,t]/self.data.COP[h]) for h in heat_pump for t in time) )

#        # enhanced MP

        P_0_max = {(h,t):max([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans]) for h in CHP for t in time}
        P_0_min = {(h,t):min([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans]) for h in CHP for t in time}
        P_0_average = {(h,t):sum([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP for t in time}
                    
        P_max = {(h,t):max([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans]) for h in CHP_sorted['ex'] for t in time}
        P_min = {(h,t):min([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans]) for h in CHP_sorted['ex']  for t in time}
        P_average = {(h,t):sum([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP_sorted['ex']  for t in time}

        Q_max = {(h,t):max([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans]) for h in CHP+heat_only for t in time}
        Q_min = {(h,t):min([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans]) for h in CHP+heat_only for t in time}
        Q_average = {(h,t):sum([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP+heat_only for t in time}


        HP_load_max = {(h,t):max([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans]) for h in heat_pump for t in time}
        HP_load_min = {(h,t):min([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans]) for h in heat_pump  for t in time}
        HP_load_average = {(h,t):sum([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in heat_pump for t in time}

        
        for h in CHP_sorted['ex']:
            for t in time:


                self.constraints.ratio_MP[h,t].rhs = P_average[h,t] + P_0_average[h,t] - self.data.r_min[h]*Q_average[h,t]
         
                self.constraints.maxprod_MP[h,t].rhs = self.data.CHP_maxprod[h] - (P_average[h,t]+ P_0_average[h,t])*self.data.rho_elec[h] - self.data.rho_heat[h]*Q_average[h,t]
                self.constraints.minprod_MP[h,t].rhs =  - Q_average[h,t]
                

                                    
        for h in heat_pump:
            for t in time:
         
                self.constraints.maxprod_MP[h,t].rhs = HP_load_average[h,t]  


#        #solve MP and update lb
        self.model.optimize()        
        self._update_lb()
        self.data.ub = self.data.z_sub + self.data.z_master
        
#        
#        # Update fixed variables for submodels 1 and solve.
        [self.submodels1[s].update_fixed_vars() for s in self.data.scenario_kmeans]
        [self.submodels1[s].optimize() for s in self.data.scenario_kmeans]
           
        [self.submodels2[s].update_fixed_vars() for s in self.data.scenario_kmeans]
        [self.submodels2[s].optimize() for s in self.data.scenario_kmeans]
#        
#        # Update bounds based on submodel rebuild
        self._update_ub()        
        self._save_vars()
#        
#        self._do_benders_step()
#        
#        # Build cuts until we reach absolute and relative tolerance,
#        # or max_iters cuts have been generated.
#        
        while (
            (np.absolute(self.data.ub - self.data.lb)>self.data.epsilon and
                len(self.data.cutlist) < self.data.max_iters)):
            # Generate new cut.
            self._do_benders_step()       


        
    def _load_data(self,W_max):
        
        #indexes
        self.data.S_kmeans = S_kmeans
        self.data.S_all_kmeans = S_all_kmeans
        self.data.scenario_load_kmeans=scenario_load_kmeans
        self.data.scenario_wind_kmeans=scenario_wind_kmeans
        self.data.scenario_supply_kmeans=scenario_supply_kmeans
        self.data.scenario_kmeans=scenario_kmeans
        self.data.scenario_list_kmeans=scenario_list_kmeans
        self.data.scenario_dic_kmeans=scenario_dic_kmeans
    
        self.data.time = time
        self.data.gen=gen
        self.data.heat_storage=heat_storage
        self.data.elec_storage=elec_storage
        self.data.heat_pump =heat_pump
        self.data.wind=wind
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP  

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

        # LOADS
        self.data.heat_load = heat_load
        
        # Elec station parameters
        self.data.elec_maxprod = elec_maxprod
        self.data.W_max = W_max
        #self.data.big_M = 1000

        self.data.elec_load = {(s,t):elec_load_scenario_kmeans[self.data.scenario_dic_kmeans[s][0],t] for s in self.data.scenario_kmeans for t in time}  
        self.data.wind_scenario = {(s,w,t):wind_scenario_kmeans[self.data.scenario_dic_kmeans[s][1],w,t] for s in self.data.scenario_kmeans for w in wind for t in time} #SCENARIOS REALIZATIONS!!!!!!!!!!!!!!!!!
        self.data.alpha = {(s,g,t):alpha_scenario_kmeans[self.data.scenario_dic_kmeans[s][2],g,t] for g in CHP+heat_only+gen+wind for s in self.data.scenario_kmeans for t in time} #BIDS
        self.data.alpha_up = {(s,g,t):alpha_up_kmeans[self.data.scenario_dic_kmeans[s][2],g,t] for g in CHP+heat_only for s in self.data.scenario_kmeans for t in time}
        self.data.alpha_down = {(s,g,t):alpha_down_kmeans[self.data.scenario_dic_kmeans[s][2],g,t] for g in CHP+heat_only for s in self.data.scenario_kmeans for t in time}       
        self.data.alpha_min = -300
        self.data.alpha_HP = 100

        #warm start
        self.data.storage_discharge_hierarchical_0 = {(h,t):storage_discharge_hierarchical_0[W_max,h,t] for h in heat_storage for t in time} 
        self.data.storage_charge_hierarchical_0 = {(h,t):storage_charge_hierarchical_0[W_max,h,t] for h in heat_storage for t in time} 
        self.data.storage_energy_hierarchical_0 = {(h,t):storage_energy_hierarchical_0[W_max,h,t] for h in heat_storage for t in time}  
        self.data.Q_hierarchical_0 = {(h,t):Q_hierarchical_0[W_max,h,t] for h in heat_pump+CHP+heat_only for t in time}               
        self.data.HP_max_load_hierarchical_0 = {(h,t):Q_hierarchical_0[W_max,h,t]/COP[h] for h in heat_pump for t in time}
        
    def _init_benders_params(self, epsilon, max_iters):

        self.data.cut_rhs = {}
        self.data.betas = []
        self.data.Qs = []
        self.data.storage_charges = []             
        self.data.storage_discharges = []
        self.data.storage_energys = []
        self.data.HP_max_loads = []
        self.data.upper_bounds = []
        self.data.lower_bounds = []
        self.data.solvetime = []        
        self.data.cutlist  =[]
        self.data.ub = gb.GRB.INFINITY
        self.data.lb = -gb.GRB.INFINITY
        self.data.max_iters = max_iters
        self.data.epsilon = epsilon
        
    def _build_model(self):
        
        self.model = gb.Model()
        self._build_variables()
        #self._build_objective()
        self._build_constraints()
        self.model.update()

    def _do_benders_step(self):
        scenario_kmeans = self.data.scenario_kmeans
        #add cut tot he master problem
        self._add_cut() #m.update() included in self._add_cut()
        #solve master pb
        self.model.optimize()
        self._update_lb() #need to m.update()?? NO
        #solve subproblems 1 & 2
        [self.submodels1[s].update_fixed_vars() for s in scenario_kmeans]
        [self.submodels1[s].optimize() for s in scenario_kmeans]
        
        [self.submodels2[s].update_fixed_vars() for s in scenario_kmeans]
        [self.submodels2[s].optimize() for s in scenario_kmeans]
        # Update bounds based on submodel rebuild
        if (self.data.z_master + sum([self.submodels1[s].model.ObjVal for s in scenario_kmeans])/S_all_kmeans < self.data.ub - 0.001*(self.data.ub - self.data.lb)):
            self._update_ub() #need to m.update()?? NO
        #self._update_ub()
        self._save_vars() #need to m.update()?? NO
        #self.model.update()
    
    def _build_variables(self):

        #indexes shortcuts 
        time = self.data.time
        #S = self.data.S
        heat_storage=self.data.heat_storage
        #elec_storage=self.data.elec_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        CHP=self.data.CHP   
        m = self.model
        
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
                      
        self.variables.beta = m.addVar(lb=-100000)
        
        m.update()
            
    def _build_constraints(self):
        
        time = self.data.time
        CHP = self.data.CHP
        heat_only = self.data.heat_only
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        m = self.model
 
        # Upper Level Constraints
        # 1st stage constraints
        # heat balance

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
                    self.variables.storage_energy[h,t1]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,t2]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,t2]-self.data.storage_loss[h])
        

            self.constraints.storage_init[h]=m.addConstr(
                self.variables.storage_energy[h,time[0]],
                gb.GRB.EQUAL,
                self.data.storage_energy_init[h]-self.data.storage_discharge_eff[h]*self.variables.storage_discharge[h,time[0]]+self.data.storage_charge_eff[h]*self.variables.storage_charge[h,time[0]]-self.data.storage_loss[h])


            self.constraints.storage_final[h]=m.addConstr(
                self.variables.storage_energy[h,time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.data.storage_energy_init[h])

        
        self.constraints.cuts = {}

#        #enhanced MP                  
                    
        self.constraints.ratio_MP = {}
        
        for h in CHP_sorted['ex']:
            for t in time:
         
                self.constraints.ratio_MP[h,t] = self.model.addConstr(
                    self.data.r_min[h]*self.variables.Q[h,t],
                    gb.GRB.LESS_EQUAL,
                    self.data.r_min[h]*self.data.heat_maxprod[h])

        self.constraints.maxprod_MP = {}
        
        for h in CHP_sorted['ex']:
            for t in time:
         
                self.constraints.maxprod_MP[h,t] = self.model.addConstr(
                    self.data.rho_heat[h]*self.variables.Q[h,t],
                    gb.GRB.LESS_EQUAL,
                    self.data.rho_heat[h]*self.data.heat_maxprod[h])

        self.constraints.minprod_MP = {}
        
        for h in CHP_sorted['ex']:
            for t in time:
         
                self.constraints.minprod_MP[h,t] = self.model.addConstr(
                    self.variables.Q[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)
                                      
        for h in heat_pump:
            for t in time:
         
                self.constraints.maxprod_MP[h,t] = self.model.addConstr(
                    self.variables.HP_max_load[h,t],
                    gb.GRB.GREATER_EQUAL,
                    0)    

                    
        m.update()
                
        
    def _add_cut(self):


        time = self.data.time
        CHP = self.data.CHP
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        m = self.model

        scenario_kmeans=self.data.scenario_kmeans
        S_all_kmeans = self.data.S_all_kmeans


        cut = len(self.data.cutlist)
        self.data.cutlist.append(cut)

        # Get sensitivities from subproblems 2
        self.data.sens_Q = {(h,t): sum([self.submodels2[s].constraints.fix_Q[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in CHP+heat_only for t in time}
        self.data.sens_storage_discharge = {(h,t): sum([self.submodels2[s].constraints.fix_storage_discharge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_charge = {(h,t): sum([self.submodels2[s].constraints.fix_storage_charge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_energy = {(h,t): sum([self.submodels2[s].constraints.fix_storage_energy[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_HP_max_load  =  {(h,t): sum([self.submodels2[s].constraints.fix_HP_max_load[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_pump for t in time}

        # Get subproblem objectives)
        self.data.z_sub = sum([self.submodels1[s].model.ObjVal for s in scenario_kmeans])/S_all_kmeans

        # Generate cut

        self.constraints.cuts[cut] = self.model.addConstr(
            self.variables.beta,
            gb.GRB.GREATER_EQUAL,
            self.data.z_sub + gb.quicksum(self.data.sens_Q[h,t]*(self.variables.Q[h,t] - self.variables.Q[h,t].x ) for h in CHP+heat_only for t in time) + gb.quicksum(self.data.sens_storage_discharge[h,t]*(self.variables.storage_discharge[h,t] - self.variables.storage_discharge[h,t].x ) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_storage_charge[h,t]*(self.variables.storage_charge[h,t]-self.variables.storage_charge[h,t].x) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_storage_energy[h,t]*(self.variables.storage_energy[h,t]-self.variables.storage_energy[h,t].x) for h in heat_storage for t in time) + gb.quicksum(self.data.sens_HP_max_load[h,t]*(self.variables.HP_max_load[h,t]-self.variables.HP_max_load[h,t].x) for h in heat_pump for t in time) )


        P_0_max = {(h,t):max([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans]) for h in CHP for t in time}
        P_0_min = {(h,t):min([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans]) for h in CHP for t in time}
        P_0_average = {(h,t):sum([self.submodels1[s].variables.P_0[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP for t in time}
                    
        P_max = {(h,t):max([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans]) for h in CHP_sorted['ex'] for t in time}
        P_min = {(h,t):min([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans]) for h in CHP_sorted['ex']  for t in time}
        P_average = {(h,t):sum([self.submodels1[s].variables.P[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP_sorted['ex']  for t in time}

        Q_max = {(h,t):max([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans]) for h in CHP+heat_only for t in time}
        Q_min = {(h,t):min([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans]) for h in CHP+heat_only for t in time}
        Q_average = {(h,t):sum([self.submodels1[s].variables.Q_up[h,t].x - self.submodels1[s].variables.Q_down[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in CHP+heat_only for t in time}


        HP_load_max = {(h,t):max([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans]) for h in heat_pump for t in time}
        HP_load_min = {(h,t):min([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans]) for h in heat_pump  for t in time}
        HP_load_average = {(h,t):sum([self.submodels1[s].variables.HP_load[h,t].x for s in scenario_kmeans])/S_all_kmeans for h in heat_pump for t in time}



        for h in CHP_sorted['ex']:
            for t in time:


                self.constraints.ratio_MP[h,t].rhs = P_average[h,t] + P_0_average[h,t] - self.data.r_min[h]*Q_average[h,t]
         
                self.constraints.maxprod_MP[h,t].rhs = self.data.CHP_maxprod[h] - (P_average[h,t]+ P_0_average[h,t])*self.data.rho_elec[h] - self.data.rho_heat[h]*Q_average[h,t]
                self.constraints.minprod_MP[h,t].rhs =  - Q_average[h,t]
                

                                    
        for h in heat_pump:
            for t in time:
         
                self.constraints.maxprod_MP[h,t].rhs = HP_load_average[h,t] 
                
                
                     
        m.update()    

    ###
    # Update upper and lower bounds for Benders' iterations
    ###
        
    def _update_lb(self):
        
        scenario_kmeans=self.data.scenario_kmeans
        S_all_kmeans = self.data.S_all_kmeans
        
        self.data.z_master_reg = self.model.ObjVal - self.variables.beta.x
        self.data.z_master = self.data.z_master_reg -  1/2*self.data.level*(sum((self.variables.Q[h,t].x-self.data.Q_ref[h,t])*(self.variables.Q[h,t].x-self.data.Q_ref[h,t]) for h in CHP+heat_only for t in time) + sum((self.variables.HP_max_load[h,t].x-self.data.HP_max_load_ref[h,t])*(self.variables.HP_max_load[h,t].x-self.data.HP_max_load_ref[h,t]) for h in heat_pump for t in time) + sum((self.variables.storage_energy[h,t].x-self.data.storage_energy_ref[h,t])*(self.variables.storage_energy[h,t].x-self.data.storage_energy_ref[h,t]) + (self.variables.storage_charge[h,t].x-self.data.storage_charge_ref[h,t])*(self.variables.storage_charge[h,t].x-self.data.storage_charge_ref[h,t]) + (self.variables.storage_discharge[h,t].x-self.data.storage_discharge_ref[h,t])*(self.variables.storage_discharge[h,t].x-self.data.storage_discharge_ref[h,t]) for h in heat_storage for t in time) )
        self.data.lb = self.data.z_master + self.variables.beta.x

        self.data.lower_bounds.append(self.data.lb)
            
    def _update_ub(self):
        
        time = self.data.time
        CHP = self.data.CHP
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        scenario_kmeans=self.data.scenario_kmeans
        S_all_kmeans = self.data.S_all_kmeans
        m = self.model 
        
        self.data.z_master_reg = self.model.ObjVal - self.variables.beta.x
        self.data.z_master = self.data.z_master_reg -  1/2*self.data.level*(sum((self.variables.Q[h,t].x-self.data.Q_ref[h,t])*(self.variables.Q[h,t].x-self.data.Q_ref[h,t]) for h in CHP+heat_only for t in time) + sum((self.variables.HP_max_load[h,t].x-self.data.HP_max_load_ref[h,t])*(self.variables.HP_max_load[h,t].x-self.data.HP_max_load_ref[h,t]) for h in heat_pump for t in time) + sum((self.variables.storage_energy[h,t].x-self.data.storage_energy_ref[h,t])*(self.variables.storage_energy[h,t].x-self.data.storage_energy_ref[h,t]) + (self.variables.storage_charge[h,t].x-self.data.storage_charge_ref[h,t])*(self.variables.storage_charge[h,t].x-self.data.storage_charge_ref[h,t]) + (self.variables.storage_discharge[h,t].x-self.data.storage_discharge_ref[h,t])*(self.variables.storage_discharge[h,t].x-self.data.storage_discharge_ref[h,t]) for h in heat_storage for t in time))
        self.data.z_sub = sum([self.submodels1[s].model.ObjVal for s in scenario_kmeans])/S_all_kmeans
        
        self.data.ub =  self.data.z_master + self.data.z_sub


        self.data.D_sens_Q_ref = {(h,t): sum([self.submodels2[s].constraints.fix_Q[h,t].pi for s in scenario_kmeans])/S_all_kmeans - self.data.sens_Q_ref[h,t] for h in CHP+heat_only for t in time}
        self.data.D_sens_storage_discharge_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_discharge[h,t].pi for s in scenario_kmeans])/S_all_kmeans - self.data.sens_storage_discharge_ref[h,t] for h in heat_storage for t in time}
        self.data.D_sens_storage_charge_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_charge[h,t].pi for s in scenario_kmeans])/S_all_kmeans - self.data.sens_storage_charge_ref[h,t] for h in heat_storage for t in time}
        self.data.D_sens_storage_energy_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_energy[h,t].pi for s in scenario_kmeans])/S_all_kmeans - self.data.sens_storage_energy_ref[h,t] for h in heat_storage for t in time}
        self.data.D_sens_HP_max_load_ref  =  {(h,t): sum([self.submodels2[s].constraints.fix_HP_max_load[h,t].pi for s in scenario_kmeans])/S_all_kmeans - self.data.sens_HP_max_load_ref[h,t] for h in heat_pump for t in time}

        self.data.sens_Q_ref = {(h,t): sum([self.submodels2[s].constraints.fix_Q[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in CHP+heat_only for t in time}
        self.data.sens_storage_discharge_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_discharge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_charge_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_charge[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_storage_energy_ref = {(h,t): sum([self.submodels2[s].constraints.fix_storage_energy[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_storage for t in time}
        self.data.sens_HP_max_load_ref  =  {(h,t): sum([self.submodels2[s].constraints.fix_HP_max_load[h,t].pi for s in scenario_kmeans])/S_all_kmeans for h in heat_pump for t in time}

        self.data.D_Q_ref = {(h,t) : self.variables.Q[h,t].x - self.data.Q_ref[h,t]  for h in CHP+heat_only for t in time}
        self.data.D_HP_max_load_ref = {(h,t) : self.variables.HP_max_load[h,t].x - self.data.HP_max_load_ref[h,t]  for h in heat_pump for t in time}
        self.data.D_storage_discharge_ref = {(h,t) : self.variables.storage_discharge[h,t].x - self.data.storage_discharge_ref[h,t] for h in heat_storage for t in time}
        self.data.D_storage_charge_ref = {(h,t) : self.variables.storage_charge[h,t].x - self.data.storage_charge_ref[h,t] for h in heat_storage for t in time}
        self.data.D_storage_energy_ref = {(h,t) : self.variables.storage_energy[h,t].x - self.data.storage_energy_ref[h,t] for h in heat_storage for t in time}
        
        self.data.Q_ref = {(h,t) : self.variables.Q[h,t].x for h in CHP+heat_only for t in time}
        self.data.HP_max_load_ref = {(h,t) : self.variables.HP_max_load[h,t].x for h in heat_pump for t in time}
        self.data.storage_discharge_ref = {(h,t) : self.variables.storage_discharge[h,t].x for h in heat_storage for t in time}
        self.data.storage_charge_ref = {(h,t) : self.variables.storage_charge[h,t].x for h in heat_storage for t in time}
        self.data.storage_energy_ref = {(h,t) : self.variables.storage_energy[h,t].x for h in heat_storage for t in time}

        self.data.level = min(1000,1/max(0.00001,(1/max(0.00001,self.data.level) + max(0,(sum([self.data.D_sens_Q_ref[h,t]*self.data.D_Q_ref[h,t] for h in CHP+heat_only for t in time]) + sum([self.data.D_sens_HP_max_load_ref[h,t]*self.data.D_HP_max_load_ref[h,t] for h in heat_pump for t in time]) + sum([self.data.D_sens_storage_discharge_ref[h,t]*self.data.D_storage_discharge_ref[h,t] for h in heat_storage for t in time]) + sum([self.data.D_sens_storage_charge_ref[h,t]*self.data.D_storage_charge_ref[h,t] for h in heat_storage for t in time]) + sum([self.data.D_sens_storage_energy_ref[h,t]*self.data.D_storage_energy_ref[h,t] for h in heat_storage for t in time]))/max(0.00001,(sum([self.data.D_sens_Q_ref[h,t]*self.data.D_sens_Q_ref[h,t] for h in CHP+heat_only for t in time]) + sum([self.data.D_sens_HP_max_load_ref[h,t]*self.data.D_sens_HP_max_load_ref[h,t] for h in heat_pump for t in time]) + sum([self.data.D_sens_storage_discharge_ref[h,t]*self.data.D_sens_storage_discharge_ref[h,t] for h in heat_storage for t in time]) + sum([self.data.D_sens_storage_charge_ref[h,t]*self.data.D_sens_storage_charge_ref[h,t] for h in heat_storage for t in time]) + sum([self.data.D_sens_storage_energy_ref[h,t]*self.data.D_sens_storage_energy_ref[h,t] for h in heat_storage for t in time])))))))

        m.setObjective(gb.quicksum(self.data.alpha[s,h,t]*self.data.rho_heat[h]*self.variables.Q[h,t] for s in scenario_kmeans for h in CHP+heat_only for t in time)/S_all_kmeans + self.variables.beta + 1/2*self.data.level*(gb.quicksum((self.variables.Q[h,t]-self.data.Q_ref[h,t])*(self.variables.Q[h,t]-self.data.Q_ref[h,t]) for h in CHP+heat_only for t in time) + gb.quicksum((self.variables.HP_max_load[h,t]-self.data.HP_max_load_ref[h,t])*(self.variables.HP_max_load[h,t]-self.data.HP_max_load_ref[h,t]) for h in heat_pump for t in time) + gb.quicksum((self.variables.storage_energy[h,t]-self.data.storage_energy_ref[h,t])*(self.variables.storage_energy[h,t]-self.data.storage_energy_ref[h,t]) + (self.variables.storage_charge[h,t]-self.data.storage_charge_ref[h,t])*(self.variables.storage_charge[h,t]-self.data.storage_charge_ref[h,t]) + (self.variables.storage_discharge[h,t]-self.data.storage_discharge_ref[h,t])*(self.variables.storage_discharge[h,t]-self.data.storage_discharge_ref[h,t]) for h in heat_storage for t in time) ),   
            gb.GRB.MINIMIZE)

        self.data.upper_bounds.append(self.data.ub)

        
    def _save_vars(self):

        time = self.data.time
        CHP = self.data.CHP
        heat_storage=self.data.heat_storage
        heat_pump=self.data.heat_pump
        heat_only=self.data.heat_only
        
        self.data.betas.append(self.variables.beta.x)
        self.data.Qs.append([self.variables.Q[h,t].x for h in CHP+heat_only for t in time])
        self.data.storage_discharges.append([self.variables.storage_discharge[h,t].x for h in heat_storage for t in time]) 
        self.data.storage_charges.append([self.variables.storage_charge[h,t].x for h in heat_storage for t in time]) 
        self.data.storage_energys.append([self.variables.storage_energy[h,t].x for h in heat_storage for t in time])                
        self.data.HP_max_loads.append([self.variables.HP_max_load[h,t].x for h in heat_pump for t in time])   
        
#%% SOLVE THE MPEC                         

#Q_hierarchical = {}
#elec_maxprod_hierarchical = {}
#elec_minprod_hierarchical = {}
#heat_cost_0_hierarchical = {}


for W_max in W_range:                 
    heat_dispatch_BD = Benders_master(W_max) 
    heat_dispatch_BD.model.params.OutputFlag = 0
#    tstart = ttt.time()
    heat_dispatch_BD.master_optimize()
#    tend = ttt.time()
#    exectime_BD[S_kmeans]=tend-tstart
#    iterations_BD[W_max] = len(heat_dispatch_BD.data.lower_bounds)
#    lb_BD[W_max] = heat_dispatch_BD.data.lb
#    ub_BD[W_max] = heat_dispatch_BD.data.ub
 
for W_max in W_range:  
    
    for t in time:
        
        for g in CHP+heat_only:
            Q_hierarchical[W_max,g,t]=heat_dispatch_BD.variables.Q[g,t].x
            
        for g in heat_pump:
            Q_hierarchical[W_max,g,t]=heat_dispatch_BD.variables.HP_max_load[g,t].x*COP[g]

        for g in heat_storage:
            Q_hierarchical[W_max,g,t]=heat_dispatch_BD.variables.storage_discharge[g,t].x - heat_dispatch_BD.variables.storage_charge[g,t].x   

    for t in time:
        
        for h in CHP_sorted['bp']:
            elec_maxprod_hierarchical[W_max,h,t] = 0
    
        for h in CHP_sorted['ex']:
            elec_maxprod_hierarchical[W_max,h,t] = (CHP_maxprod[h]-rho_heat[h]*Q_hierarchical[W_max,h,t])/rho_elec[h] - r_min[h]*Q_hierarchical[W_max,h,t]    

        for g in gen:
            elec_maxprod_hierarchical[W_max,g,t] = elec_maxprod[g]   

        for g in wind:
            elec_maxprod_hierarchical[W_max,g,t] = W_max 
            
        for g in heat_pump:
            elec_maxprod_hierarchical[W_max,g,t] = Q_hierarchical[W_max,g,t]/COP[g]
            
        for h in CHP:
            elec_minprod_hierarchical[W_max,h,t] = r_min[h]*Q_hierarchical[W_max,h,t]


    heat_cost_0_hierarchical[W_max]=heat_dispatch_BD.data.ub
       
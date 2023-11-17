from copy import deepcopy
import HX_module as HX
import COMPAND_module as CP
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as PLT
import numpy as np
# test 용 import
from HP_dataclass import ProcessFluid, Settings, Outputs


class VCHP():
    def __init__(self, InCond, OutCond, InEvap, OutEvap, inputs):
        self.InCond = InCond
        self.OutCond = OutCond
        self.InEvap = InEvap
        self.OutEvap = OutEvap
        self.inputs = inputs
    
    def __call__(self):
        outputs = Outputs()
        
        InCond_REF = ProcessFluid(Y=self.inputs.Y)
        OutCond_REF = ProcessFluid(Y=self.inputs.Y)
        InEvap_REF = ProcessFluid(Y=self.inputs.Y)
        OutEvap_REF = ProcessFluid(Y=self.inputs.Y)
        
        p_crit = PropsSI('PCRIT','',0,'',0, list(self.inputs.Y.keys())[0])
        T_crit = PropsSI('TCRIT','',0,'',0, list(self.inputs.Y.keys())[0])
        
        InCond_REF.p_crit = p_crit
        OutCond_REF.p_crit = p_crit
        InEvap_REF.p_crit = p_crit
        OutEvap_REF.p_crit = p_crit
        
        InCond_REF.T_crit = T_crit
        OutCond_REF.T_crit = T_crit
        InEvap_REF.T_crit = T_crit
        OutEvap_REF.T_crit = T_crit
        
        (self.InCond, self.OutCond, self.InEvap, self.OutEvap, no_input) = self.Input_Processing(self.InCond, self.OutCond, self.InEvap, self.OutEvap, self.inputs)
        evap_ph = 0
        cond_ph = 0
        
        if inputs.layout == 'bas':
            (self.InCond, self.OutCond, self.InEvap, self.OutEvap, self.InCond_REF, self.OutCond_REF, self.InEvap_REF, self.OutEvap_REF, outputs) = self.Cycle_Solver(self.InCond, self.OutCond, self.InEvap, self.OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, self.inputs, outputs, no_input, cond_ph, evap_ph)
        
        if inputs.layout == '3eco':
            (self.InCond, self.OutCond, self.InEvap, self.OutEvap, self.InCond_REF, self.OutCond_REF, self.InEvap_REF, self.OutEvap_REF, outputs) = self.Three_Economizer_Solver(self.InCond, self.OutCond, self.InEvap, self.OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, self.inputs, outputs, no_input, cond_ph, evap_ph)
        
        #self.Post_Processing(outputs)
        #self.Plot_diagram(self.InCond_REF, self.OutCond_REF, self.InEvap_REF, self.OutEvap_REF, self.inputs, outputs)
        
        return (self.InCond, self.OutCond, self.InEvap, self.OutEvap, self.InCond_REF, self.OutCond_REF, self.InEvap_REF, self.OutEvap_REF, outputs)
        
    def Input_Processing(self, InCond, OutCond, InEvap, OutEvap, inputs):
        no_InCondT = 0
        no_Condm = 0
        no_OutCondT = 0
        no_InEvapT = 0
        no_Evapm = 0
        no_OutEvapT = 0
        
        
        if inputs.second == 'steam':
            (InCond, OutCond) = self.Steam_module(InCond, OutCond, inputs)
        
        else: # 일반 공정
            if InCond.p <= 0.0:
                InCond.p = 101300.0
        
            if OutCond.p <= 0.0:
             OutCond.p = 101300.0
             
            if InCond.T <= 0.0:
                no_InCondT = 1
        
            if InCond.m <= 0 and OutCond.m <= 0 :
                no_Condm = 1
            else:
                if InCond.m == 0:
                    InCond.m = OutCond.m # shallow copy
                else:
                    OutCond.m = InCond.m 
                
            if OutCond.T <= 0.0:
                no_OutCondT = 1
            
            
        if InEvap.p <= 0.0:
            InEvap.p = 101300.0
        
        if OutEvap.p <= 0.0:
            OutEvap.p = 101300.0
        
        
        if InEvap.T <= 0.0:
            no_InEvapT = 1
        
        if InEvap.m <= 0 and OutEvap.m <= 0 :
            no_Evapm = 1
        else:
            if InEvap.m == 0:
                InEvap.m = OutEvap.m
            else:
                OutEvap.m = InEvap.m
                
        if OutEvap.T <= 0.0:
            no_OutEvapT = 1
            
        no_inputs_sum = no_InCondT + no_Condm + no_OutCondT + no_InEvapT + no_Evapm + no_OutEvapT
        
        if no_inputs_sum == 0:
            no_input = 'Overdefine'
            return (InCond, OutCond, InEvap, OutEvap, no_input)
        elif no_inputs_sum > 1:
            no_input = 'Underdefine'
            return (InCond, OutCond, InEvap, OutEvap, no_input)
        else:    
            if no_InCondT == 1:
                OutCond.h = PropsSI('H','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutCond.Cp = PropsSI('C','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                InEvap.h = PropsSI('H','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InEvap.Cp = PropsSI('C','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                OutEvap.h = PropsSI('H','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                InEvap.q = (OutEvap.h - InEvap.h)*InEvap.m
                OutEvap.q = InEvap.q 
                no_input = 'InCondT'
            elif no_OutCondT == 1:
                InCond.h = PropsSI('H','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InCond.Cp = PropsSI('C','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InEvap.h = PropsSI('H','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InEvap.Cp = PropsSI('C','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                OutEvap.h = PropsSI('H','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                InEvap.q = (OutEvap.h - InEvap.h)*InEvap.m
                OutEvap.q = InEvap.q 
                no_input = 'OutCondT'
            elif no_Condm == 1:
                InCond.h = PropsSI('H','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InCond.Cp = PropsSI('C','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                OutCond.h = PropsSI('H','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutCond.Cp = PropsSI('C','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                InEvap.h = PropsSI('H','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InEvap.Cp = PropsSI('C','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                OutEvap.h = PropsSI('H','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                InEvap.q = (OutEvap.h - InEvap.h)*InEvap.m
                OutEvap.q = InEvap.q 
                no_input = 'Condm'
            
            elif no_InEvapT == 1:
                InCond.h = PropsSI('H','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InCond.Cp = PropsSI('C','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                OutCond.h = PropsSI('H','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutCond.Cp = PropsSI('C','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutEvap.h = PropsSI('H','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                InCond.q = (OutCond.h - InCond.h)*InCond.m
                OutCond.q = InCond.q
                no_input = 'InEvapT'
            elif no_OutEvapT == 1:
                InCond.h = PropsSI('H','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InCond.Cp = PropsSI('C','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                OutCond.h = PropsSI('H','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutCond.Cp = PropsSI('C','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                InEvap.h = PropsSI('H','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InEvap.Cp = PropsSI('C','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InCond.q = (OutCond.h - InCond.h)*InCond.m
                OutCond.q = InCond.q
                no_input = 'OutEvapT'
            elif no_Evapm == 1:
                InCond.h = PropsSI('H','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                InCond.Cp = PropsSI('C','T',InCond.T, 'P',InCond.p, InCond.fluidmixture)
                OutCond.h = PropsSI('H','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                OutCond.Cp = PropsSI('C','T',OutCond.T, 'P',OutCond.p, OutCond.fluidmixture)
                InEvap.h = PropsSI('H','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                InEvap.Cp = PropsSI('C','T',InEvap.T, 'P',InEvap.p, InEvap.fluidmixture)
                OutEvap.h = PropsSI('H','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P',OutEvap.p, OutEvap.fluidmixture)
                InCond.q = (OutCond.h - InCond.h)*InCond.m
                OutCond.q = InCond.q
                no_input = 'Evapm'
                
        return (InCond, OutCond, InEvap, OutEvap, no_input)
    
    def Steam_module(self, InCond, OutCond, inputs):
        p_flash = PropsSI('P','T',inputs.T_steam,'Q',1.0, InCond.fluidmixture)
        OutCond.p = PropsSI('P','T',OutCond.T+0.0001, 'Q', 0.0, InCond.fluidmixture)
        OutCond.h = PropsSI('H','T',OutCond.T+0.0001, 'Q', 0.0, InCond.fluidmixture)
        X_flash = PropsSI('Q','H',OutCond.h,'P',p_flash, InCond.fluidmixture)
        OutCond.m = inputs.m_steam / X_flash
        InCond.m = OutCond.m
        m_sat_liq = (1-X_flash)*OutCond.m
        h_sat_liq = PropsSI('H','P',p_flash,'Q',0.0, InCond.fluidmixture)
        h_makeup = PropsSI('H','T',inputs.T_makeup,'P',p_flash, InCond.fluidmixture)
        
        InCond.h = (m_sat_liq*h_sat_liq + inputs.m_makeup*h_makeup)/OutCond.m
        InCond.T = PropsSI('T','H',InCond.h,'P',p_flash, InCond.fluidmixture)
        InCond.p = OutCond.p
        return (InCond, OutCond)
        
    def HighPressure_Solver(self, InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph):
        if inputs.cycle == 'scc':
            cond_p_ub = min(5*InCond_REF.p_crit, 1.0e8)
            cond_p_lb = InCond_REF.p_crit
        else:
            cond_p_ub = InCond_REF.p_crit
            if no_input == 'InCondT':
                cond_p_lb = PropsSI('P','T',OutCond.T,'Q',1.0,OutCond_REF.fluidmixture)
            else:
                cond_p_lb = PropsSI('P','T',InCond.T,'Q',1.0,OutCond_REF.fluidmixture)

        cond_a = 1
        while cond_a:
            InCond_REF.p = 0.5*(cond_p_ub+cond_p_lb)
            OutCond_REF.p = InCond_REF.p*(1-inputs.cond_dp)
            
            
            OutCond_REF.T = PropsSI('T','P',OutCond_REF.p,'Q',0.0,OutCond_REF.fluidmixture) - inputs.DSC
            if inputs.DSC == 0:
                OutCond_REF.h = PropsSI('H','P',OutCond_REF.p,'Q',0.0,OutCond_REF.fluidmixture)
                if inputs.expand_eff > 0.0:
                    OutCond_REF.s = PropsSI('S','P',OutCond_REF.p,'Q',0.0,OutCond_REF.fluidmixture)
            else:
                OutCond_REF.h = PropsSI('H','T',OutCond_REF.T,'P',OutCond_REF.p, OutCond_REF.fluidmixture)
                if inputs.expand_eff > 0.0:
                    OutCond_REF.s = PropsSI('S','T',OutCond_REF.T,'P',OutCond_REF.p, OutCond_REF.fluidmixture)
                
            if inputs.layout == 'bas':
                InComp = OutEvap_REF
                InExpand = OutCond_REF
                
                comp = CP.Compander_module(InComp, InCond_REF)
                (inputs.DSH, cond_a) = comp.COMP(eff_isen = inputs.comp_eff, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                InCond_REF = comp.primary_out
                
                expand = CP.Compander_module(InExpand, InEvap_REF)
                expand.EXPAND(eff_isen = inputs.expand_eff, eff_mech = inputs.mech_eff)
                InEvap_REF = expand.primary_out
                
            elif inputs.layout == '3eco':
                outputs.eco1_p = InEvap_REF.p*(OutCond_REF.p/InEvap_REF.p)**self.eco1_frac
                eco1_h_vap = PropsSI("H","P",outputs.eco1_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco1_h_liq = PropsSI("H","P",outputs.eco1_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                outputs.eco2_p = outputs.eco1_p*(OutCond_REF.p/InEvap_REF.p)**self.eco2_frac
                eco2_h_vap = PropsSI("H","P",outputs.eco2_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco2_h_liq = PropsSI("H","P",outputs.eco2_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                outputs.eco3_p = outputs.eco2_p*(OutCond_REF.p/InEvap_REF.p)**self.eco3_frac    
                eco3_h_vap = PropsSI("H","P",outputs.eco3_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco3_h_liq = PropsSI("H","P",outputs.eco3_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                outputs.eco3_x = (OutCond_REF.h - eco3_h_liq)/(eco3_h_vap - eco3_h_liq)
                outputs.eco2_x = (eco3_h_liq - eco2_h_liq)/(eco2_h_vap - eco2_h_liq)
                outputs.eco1_x = (eco2_h_liq - eco1_h_liq)/(eco1_h_vap - eco1_h_liq)
                                
                Out_LPcompS1 = deepcopy(OutEvap_REF)
                Out_LPcompS1.p = outputs.eco1_p
                
                LPcompS1 = CP.Compander_module(OutEvap_REF, Out_LPcompS1)
                (inputs.DSH, cond_a) = LPcompS1.COMP(eff_isen = inputs.eff_LPcompS1, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                Out_LPcompS1 = LPcompS1.primary_out
                
                In_LPcompS2 = deepcopy(Out_LPcompS1)
                In_LPcompS2.h = Out_LPcompS1.h*(1.0 - outputs.eco1_x)+eco1_h_vap*outputs.eco1_x
                In_LPcompS2.T = PropsSI("T","H",In_LPcompS2.h,"P",In_LPcompS2.p,In_LPcompS2.fluidmixture)
                In_LPcompS2.s = PropsSI("S","T",In_LPcompS2.T,"P",In_LPcompS2.p,In_LPcompS2.fluidmixture)
                
                Out_LPcompS2 = deepcopy(In_LPcompS2)
                Out_LPcompS2.p = outputs.eco2_p
                LPcompS2 = CP.Compander_module(In_LPcompS2, Out_LPcompS2)
                (dummy_dsh, dummy_a) = LPcompS2.COMP(eff_isen = inputs.eff_LPcompS2, eff_mech = inputs.mech_eff, DSH = 0.0)
                Out_LPcompS2 = LPcompS2.primary_out
                                    
                In_HPcompS1 = deepcopy(Out_LPcompS2)
                In_HPcompS1.h = Out_LPcompS2.h*(1.0 - outputs.eco2_x) + eco2_h_vap*outputs.eco2_x
                In_HPcompS1.T = PropsSI("T","H",In_HPcompS1.h,"P",In_HPcompS1.p,In_HPcompS1.fluidmixture)
                In_HPcompS1.s = PropsSI("S","T",In_HPcompS1.T,"P",In_HPcompS1.p,In_HPcompS1.fluidmixture)
                
                Out_HPcompS1 = deepcopy(In_HPcompS1)
                Out_HPcompS1.p = outputs.eco3_p
                HPcompS1 = CP.Compander_module(In_HPcompS1, Out_HPcompS1)
                (dummy_dsh, dummy_a) = HPcompS1.COMP(eff_isen = inputs.eff_HPcompS1, eff_mech = inputs.mech_eff, DSH = 0.0)
                Out_HPcompS1 = HPcompS1.primary_out
                
                In_HPcompS2 = deepcopy(Out_HPcompS1)
                In_HPcompS2.h = Out_HPcompS1.h*(1.0 - outputs.eco3_x) + eco3_h_vap*outputs.eco3_x
                In_HPcompS2.T = PropsSI("T","H",In_HPcompS2.h,"P",In_HPcompS2.p, In_HPcompS2.fluidmixture)
                In_HPcompS2.s = PropsSI("S","T",In_HPcompS2.T,"P",In_HPcompS2.p, In_HPcompS2.fluidmixture)
                
                HPcompS2 = CP.Compander_module(In_HPcompS2, InCond_REF)
                (dummy_dsh, dummy_a) = HPcompS2.COMP(eff_isen = inputs.eff_HPcompS2, eff_mech = inputs.mech_eff, DSH = 0.0)
                InCond_REF = HPcompS2.primary_out
                
                InEvap_REF.h = eco1_h_liq
                InEvap_REF.T = PropsSI("T","H",InEvap_REF.h,"P",InEvap_REF.p,InEvap_REF.fluidmixture)
                
                
            if (no_input == 'InCondT') or (no_input == 'OutCondT') or (no_input == 'Condm'):
                InEvap_REF.m = InEvap.q/(InEvap_REF.h - OutEvap_REF.h)
                OutEvap_REF.m = InEvap_REF.m
                
                if inputs.layout == 'bas':
                    InCond_REF.m = InEvap_REF.m
                    OutCond_REF.m = InEvap_REF.m
                    outputs.Wcomp = comp.Pspecific*InCond_REF.m
                    outputs.Wexpand = expand.Pspecific*InEvap_REF.m
                elif inputs.layout == '3eco':
                    InCond_REF.m = InEvap_REF.m/(1.0-outputs.eco1_x)/(1.0-outputs.eco2_x)/(1.0-outputs.eco3_x)
                    OutCond_REF.m = InCond_REF.m
                    outputs.W_LPcompS1 = LPcompS1.Pspecific*OutEvap_REF.m
                    outputs.W_LPcompS2 = LPcompS2.Pspecific*OutEvap_REF.m/(1.0-outputs.eco1_x)
                    outputs.W_HPcompS1 = HPcompS1.Pspecific*OutEvap_REF.m/(1.0-outputs.eco1_x)/(1.0-outputs.eco2_x)
                    outputs.W_HPcompS2 = HPcompS1.Pspecific*InCond_REF.m
                    
                InEvap_REF.q = -InEvap.q
                OutEvap_REF.q = -OutEvap.q
                
                InCond_REF.q = (OutCond_REF.h - InCond_REF.h)*InCond_REF.m
                OutCond_REF.q = InCond_REF.q
                InCond.q = -InCond_REF.q
                OutCond.q = -InCond_REF.q
                if no_input == 'InCondT':
                    InCond.h = OutCond.h - OutCond.q/OutCond.m
                    try:    
                        InCond.T = PropsSI('T','P',InCond.p, 'H', InCond.h, InCond.fluidmixture)
                        InCond.Cp = PropsSI('C','T',InCond.T, 'P', InCond.p, InCond.fluidmixture)
                    except:    
                        InCond.T = OutCond.T
                        while 1:
                            H_virtual = PropsSI('H','T', InCond.T, 'P', InCond.p, InCond.fluidmixture)
                            InCond.Cp = PropsSI('C','T',InCond.T, 'P', InCond.p, InCond.fluidmixture)
                            err_h = H_virtual - InCond.h
                            InCond.T = InCond.T - err_h/InCond.Cp
                            if err_h/InCond.h < 1.0e-5:
                                break
                elif no_input == 'OutCondT':
                    OutCond.h = InCond.h + InCond.q/InCond.m
                    try:    
                        OutCond.T = PropsSI('T','P',OutCond.p, 'H', OutCond.h, OutCond.fluidmixture)
                        OutCond.Cp = PropsSI('C','T',OutCond.T, 'P', OutCond.p, OutCond.fluidmixture)
                    except:
                        OutCond.T = InCond.T
                        while 1:
                            H_virtual = PropsSI('H','T', OutCond.T, 'P', OutCond.p, OutCond.fluidmixture)
                            OutCond.Cp = PropsSI('C','T',OutCond.T, 'P', OutCond.p, OutCond.fluidmixture)
                            err_h = H_virtual - OutCond.h
                            OutCond.T = OutCond.T - err_h/OutCond.Cp
                            if err_h/OutCond.h < 1.0e-5:
                                break
                elif no_input == 'Condm':
                    InCond.m = InCond.q/(OutCond.h - InCond.h)
                    OutCond.m = InCond.m
                
            elif (no_input == 'InEvapT') or (no_input == 'OutEvapT') or (no_input == 'Evapm'):
                InCond_REF.m = InCond.q/(InCond_REF.h - OutCond_REF.h)
                OutCond_REF.m = InCond_REF.m
                
                if inputs.layout == 'bas':
                    InEvap_REF.m = InCond_REF.m
                    OutEvap_REF.m = InCond_REF.m
                    outputs.Wcomp = comp.Pspecific*InCond_REF.m
                    outputs.Wexpand = expand.Pspecific*InEvap_REF.m
                elif inputs.layout == '3eco':
                    InEvap_REF.m = InCond_REF.m*(1.0-outputs.eco1_x)*(1.0-outputs.eco2_x)*(1.0-outputs.eco3_x)
                    OutEvap_REF.m = InEvap_REF.m
                    outputs.W_LPcompS1 = LPcompS1.Pspecific*OutEvap_REF.m
                    outputs.W_LPcompS2 = LPcompS2.Pspecific*OutEvap_REF.m/(1.0-outputs.eco1_x)
                    outputs.W_HPcompS1 = HPcompS1.Pspecific*OutEvap_REF.m/(1.0-outputs.eco1_x)/(1.0-outputs.eco2_x)
                    outputs.W_HPcompS2 = HPcompS1.Pspecific*InCond_REF.m
                    
                InCond_REF.q = -InCond.q
                OutCond_REF.q = -InCond.q
                
                InEvap_REF.q = (OutEvap_REF.h - InEvap_REF.h)*InEvap_REF.m
                OutEvap_REF.q = InEvap_REF.q
                InEvap.q = -InEvap_REF.q
                OutEvap.q = -InEvap_REF.q
                
                if no_input == 'InEvapT':
                    InEvap.h = OutEvap.h - OutEvap.q/OutEvap.m
                    try:    
                        InEvap.T = PropsSI('T','P',InEvap.p, 'H', InEvap.h, InEvap.fluidmixture)
                        InEvap.Cp = PropsSI('C','T',InEvap.T, 'P', InEvap.p, InEvap.fluidmixture)
                    except:    
                        InEvap.T = OutEvap.T
                        while 1:
                            H_virtual = PropsSI('H','T', InEvap.T, 'P', InEvap.p, InEvap.fluidmixture)
                            InEvap.Cp = PropsSI('C','T',InEvap.T, 'P', InEvap.p, InEvap.fluidmixture)
                            err_h = H_virtual - InEvap.h
                            InEvap.T = InEvap.T - err_h/InEvap.Cp
                            if err_h/InEvap.h < 1.0e-5:
                                break
                            
                elif no_input == 'OutEvapT':
                    OutEvap.h = InEvap.h + InEvap.q/InEvap.m
                    try:
                        OutEvap.T = PropsSI('T','P',OutEvap.p, 'H', OutEvap.h, OutEvap.fluidmixture)
                        OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P', OutEvap.p, OutEvap.fluidmixture)
                    except:
                        OutEvap.T = InEvap.T
                        while 1:
                            H_virtual = PropsSI('H','T', OutEvap.T, 'P', OutEvap.p, OutEvap.fluidmixture)
                            OutEvap.Cp = PropsSI('C','T',OutEvap.T, 'P', OutEvap.p, OutEvap.fluidmixture)
                            err_h = H_virtual - OutEvap.h
                            OutEvap.T = OutEvap.T - err_h/OutEvap.Cp
                            if err_h/OutEvap.h < 1.0e-5:
                                break
                    
                elif no_input == 'Evapm':
                    InEvap.m = InEvap.q/(OutEvap.h - InEvap.h)
                    OutEvap.m = InEvap.m
                
            cond = HX.Heatexchanger_module(InCond_REF, OutCond_REF, 1, InCond, OutCond, cond_ph)
        
            if inputs.cond_type == 'fthe':
                (outputs.cond_Tarray, outputs.cond_parray) = cond.FTHE(N_element=inputs.cond_N_element, N_row = inputs.cond_N_row)
                self.cond_err = (inputs.cond_T_lm - cond.T_lm)/inputs.cond_T_lm
                
            elif inputs.cond_type == 'phe':
                (outputs.cond_Tarray, outputs.cond_parray) = cond.PHE(N_element=inputs.cond_N_element)
                self.cond_err = (inputs.cond_T_pp - cond.T_pp)/inputs.cond_T_pp
            
            OutCond_REF = cond.primary_out
            
            if cond.T_rvs == 1:
                cond_p_lb = InCond_REF.p
            else:
                if self.cond_err < 0:
                    cond_p_ub = InCond_REF.p
                else:
                    cond_p_lb = InCond_REF.p
            
            if abs(self.cond_err) < inputs.tol:
                self.cond_conv_err = 0
                cond_a = 0
            elif (cond_p_ub - cond_p_lb)/(0.5*(cond_p_ub + cond_p_lb)) < inputs.tol:
                self.cond_conv_err = 1
                cond_a = 0
        
        outputs.cond_UA = cond.UA    
        
        return (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    
    def Cycle_Solver(self,InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph, evap_ph):
        if no_input == 'InEvapT':
            evap_p_ub = PropsSI('P','T',OutEvap.T, 'Q', 1.0, InEvap_REF.fluidmixture)        
        else:
            evap_p_ub = PropsSI('P','T',InEvap.T, 'Q', 1.0, InEvap_REF.fluidmixture)
            
        evap_p_lb = 101300.0
        evap_a = 1
        
        while evap_a: 
            OutEvap_REF.p = 0.5*(evap_p_lb+evap_p_ub)
            InEvap_REF.p = OutEvap_REF.p/(1.0-inputs.evap_dp)
            
            OutEvap_REF_Tvap = PropsSI('T','P',OutEvap_REF.p, 'Q', 1.0, OutEvap_REF.fluidmixture)
            OutEvap_REF.T = OutEvap_REF_Tvap + inputs.DSH
            if inputs.DSH == 0:
                OutEvap_REF.h = PropsSI('H','P',OutEvap_REF.p, 'Q', 1.0, OutEvap_REF.fluidmixture)
                OutEvap_REF.s = PropsSI('S','P',OutEvap_REF.p, 'Q', 1.0, OutEvap_REF.fluidmixture)
            else:
                OutEvap_REF.h = PropsSI('H','T',OutEvap_REF.T, 'P', OutEvap_REF.p ,OutEvap_REF.fluidmixture)
                OutEvap_REF.s = PropsSI('S','T',OutEvap_REF.T, 'P', OutEvap_REF.p ,OutEvap_REF.fluidmixture)
            
            (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = self.HighPressure_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph)
            
            evap = HX.Heatexchanger_module(InEvap_REF, OutEvap_REF, 1, InEvap, OutEvap, evap_ph)
            
            if inputs.evap_type == 'fthe':
                evap.FTHE(N_element = inputs.evap_N_element, N_row = inputs.evap_N_row)
                self.evap_err = (inputs.evap_T_lm - evap.T_lm)/inputs.evap_T_lm
            elif inputs.evap_type == 'phe':
                evap.PHE(N_element= inputs.evap_N_element)
                self.evap_err = (inputs.evap_T_pp - evap.T_pp)/inputs.evap_T_pp
            
            OutEvap_REF = evap.primary_out
            
            if evap.T_rvs == 1:
                evap_p_ub = OutEvap_REF.p                    
            else:
                if self.evap_err < 0:
                    evap_p_lb = OutEvap_REF.p
                else:
                    evap_p_ub = OutEvap_REF.p
                    
            if abs(self.evap_err) < inputs.tol:
                self.evap_conv_err = 0
                evap_a = 0
            elif (evap_p_ub - evap_p_lb)/(0.5*(evap_p_ub + evap_p_lb)) < inputs.tol:
                self.evap_conv_err = 1
                evap_a = 0
            
            if inputs.layout == 'bas':
                outputs.COP_heating = abs(OutCond.q)/(outputs.Wcomp - outputs.Wexpand)
            elif inputs.layout == '3eco':
                outputs.COP_heating = abs(OutCond.q)/(outputs.W_LPcompS1+outputs.W_LPcompS2+outputs.W_HPcompS1+outputs.W_HPcompS2)
                
        
        outputs.DSH = inputs.DSH
        outputs.evap_UA = evap.UA
        
        return (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    
    def Three_Economizer_Solver(self,InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph, evap_ph):
        self.eco1_frac = inputs.frac_list[1]
        self.eco2_frac = inputs.frac_list[2]
        self.eco3_frac = inputs.frac_list[3]
        
        (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = self.Cycle_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph, evap_ph)
        
        return(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    
    def Plot_diagram(self, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs):
        (p_array, h_array, T_array, s_array, p_points, h_points, s_points, T_points) = self.Dome_Draw(InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs)
        fig_ph, ax_ph = PLT.subplots()
        ax_ph.plot([i/1.0e3 for i in h_array], [i/1.0e5 for i in p_array],'b--')
        ax_ph.set_xlabel('Enthalpy [kJ/kg]',fontsize = 15)
        ax_ph.set_ylabel('Pressure [bar]',fontsize = 15)
        ax_ph.set_title('Pressure-Enthalpy Diagram\nRefrigerant:{}'.format(list(inputs.Y.keys())[0]),fontsize = 18)
        ax_ph.tick_params(axis = 'x', labelsize = 13)
        ax_ph.tick_params(axis = 'y', labelsize = 13)
        
        fig_ts, ax_ts = PLT.subplots()
        ax_ts.plot([i/1.0e3 for i in s_array], [i-273.15 for i in T_array],'b--')
        ax_ts.set_xlabel('Entropy [kJ/kg-K]',fontsize = 15)
        ax_ts.set_ylabel('Temperature [℃]',fontsize = 15)
        ax_ts.set_title('Temperature-Entropy Diagram\nRefrigerant:{}'.format(list(inputs.Y.keys())[0]),fontsize = 18)
        ax_ts.tick_params(axis = 'x', labelsize = 13)
        ax_ts.tick_params(axis = 'y', labelsize = 13)
        
        ax_ts.plot([i/1.0e3 for i in s_points], [i-273.15 for i in T_points],'bo-')
        ax_ph.plot([i/1.0e3 for i in h_points], [i/1.0e5 for i in p_points],'bo-')
        
        fig_ph.savefig('.\Figs\Ph_diagram.png',dpi=300)
        fig_ts.savefig('.\Figs\Ts_diagram.png',dpi=300)
        
    def Dome_Draw(self, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs):        
        P_crit = PropsSI('PCRIT','',0,'',0,InCond_REF.fluidmixture)
        P_trp = PropsSI('PTRIPLE','',0,'',0,InCond_REF.fluidmixture)
        T_crit = PropsSI('TCRIT','',0,'',0,InCond_REF.fluidmixture)
        H_crit = 0.5*(PropsSI('H','P',P_crit*0.999,'Q',0.0,InCond_REF.fluidmixture)+PropsSI('H','P',P_crit*0.999,'Q',1.0,InCond_REF.fluidmixture))
        S_crit = 0.5*(PropsSI('S','P',P_crit*0.999,'Q',0.0,InCond_REF.fluidmixture)+PropsSI('S','P',P_crit*0.999,'Q',1.0,InCond_REF.fluidmixture))
        
        try:
            Pliq_array = [101300.0+(P_crit*0.99 - 101300.0)*i/49 for i in range(50)]
            Pvap_array = [P_crit*0.99 - (P_crit*0.99 - 101300.0)*i/49 for i in range(50)]
            Tliq_array = [PropsSI('T','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            Tvap_array = [PropsSI('T','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
            hliq_array = [PropsSI('H','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            hvap_array = [PropsSI('H','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
            sliq_array = [PropsSI('S','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            svap_array = [PropsSI('S','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
        except:
            Pliq_array = [P_trp*1.01+(P_crit*0.99 - P_trp*1.01)*i/49 for i in range(50)]
            Pvap_array = [P_crit*0.99 - (P_crit*0.99 - P_trp*1.01)*i/49 for i in range(50)]
            Tliq_array = [PropsSI('T','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            Tvap_array = [PropsSI('T','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
            hliq_array = [PropsSI('H','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            hvap_array = [PropsSI('H','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
            sliq_array = [PropsSI('S','P',i,'Q',0.0,InCond_REF.fluidmixture) for i in Pliq_array]
            svap_array = [PropsSI('S','P',i,'Q',1.0,InCond_REF.fluidmixture) for i in Pvap_array]
            
        p_array = Pliq_array+[P_crit]+Pvap_array
        T_array = Tliq_array+[T_crit]+Tvap_array
        h_array = hliq_array+[H_crit]+hvap_array
        s_array = sliq_array+[S_crit]+svap_array
        
        OutEvap_REF_Tvap = PropsSI('T','P',OutEvap_REF.p,'Q',1.0,OutCond_REF.fluidmixture)
        OutEvap_REF_svap = PropsSI('S','P',OutEvap_REF.p,'Q',1.0,OutCond_REF.fluidmixture)
        
        
        InCond_REF_Tvap = PropsSI('T','P',InCond_REF.p,'Q',1.0,OutCond_REF.fluidmixture)
        InCond_REF_svap = PropsSI('S','P',InCond_REF.p,'Q',1.0,OutCond_REF.fluidmixture)
        OutCond_REF_Tliq = PropsSI('T','P',OutCond_REF.p,'Q',0.0,OutCond_REF.fluidmixture)
        OutCond_REF_sliq = PropsSI('S','P',OutCond_REF.p,'Q',0.0,OutCond_REF.fluidmixture)
            
        OutCond_REF.s = PropsSI('S','T',OutCond_REF.T,'P',OutCond_REF.p,OutCond_REF.fluidmixture)
        InEvap_REF.s = PropsSI('S','H',InEvap_REF.h,'P',InEvap_REF.p,InEvap_REF.fluidmixture)
        
            
        s_points = [OutEvap_REF_svap, OutEvap_REF.s, InCond_REF.s, InCond_REF_svap, OutCond_REF_sliq, OutCond_REF.s, InEvap_REF.s, OutEvap_REF_svap]
        T_points = [OutEvap_REF_Tvap, OutEvap_REF.T, InCond_REF.T, InCond_REF_Tvap, OutCond_REF_Tliq, OutCond_REF.T, InEvap_REF.T, OutEvap_REF_Tvap]
        h_points = [OutEvap_REF.h, InCond_REF.h, OutCond_REF.h, InEvap_REF.h, OutEvap_REF.h]
        p_points = [OutEvap_REF.p, InCond_REF.p, OutCond_REF.p, InEvap_REF.p, OutEvap_REF.p]
            
        return (p_array, h_array, T_array, s_array, p_points, h_points, s_points, T_points)
            
            
    def Post_Processing(self, outputs):
        print('Heating COP:{:.3f}, Cooling COP:{:.3f}'.format(outputs.COP_heating, outputs.COP_heating-1))
        print('Refrigerant:{}'.format(self.OutCond_REF.fluidmixture))
        print('Q heating: {:.3f} [kW] ({:.3f} [usRT])'.format(self.OutCond.q/1000, self.OutCond.q/3516.8525))
        print('Q cooling: {:.3f} [kW] ({:.3f} [usRT])'.format(self.OutEvap_REF.q/1000, self.OutEvap_REF.q/3516.8525))
        print('Q comp: {:.3f} [kW]'.format(outputs.Wcomp/1000))
        print('Hot fluid Inlet T:{:.3f}[℃]/P:{:.3f}[bar]/m:{:.3f}[kg/s]:   -------> Hot fluid Outlet T:{:.3f}[℃]/P:{:.3f}[bar]/m:{:.3f}[kg/s]'.format(self.InCond.T-273.15, self.InCond.p/1.0e5, self.InCond.m, self.OutCond.T-273.15, self.OutCond.p/1.0e5, self.OutCond.m))
        print('Cold fluid Outlet T:{:.3f}[℃]/P:{:.3f}[bar]/m:{:.3f}[kg/s]: <------- Cold fluid Inlet T:{:.3f}[℃]/P:{:.3f}[bar]/m:{:.3f}[kg/s]'.format(self.OutEvap.T-273.15, self.OutEvap.p/1.0e5, self.OutEvap.m, self.InEvap.T-273.15, self.InEvap.p/1.0e5, self.InEvap.m))
        print('Plow: {:.3f} [bar], Phigh: {:.3f} [bar], PR: {:.3f}[-]'.format(self.OutEvap_REF.p/1.0e5, self.InCond_REF.p/1.0e5, self.InCond_REF.p/self.OutEvap_REF.p))
        Tlow = PropsSI('T','P',0.5*(self.OutEvap_REF.p+self.InEvap_REF.p),'Q',0.5,self.OutEvap_REF.fluidmixture)
        try:
            Thigh = PropsSI('T','P',0.5*(self.OutCond_REF.p+self.InCond_REF.p),'Q',0.5,self.OutCond_REF.fluidmixture)
        except:
            Thigh = 0
        print('Tlow: {:.3f} [℃], Thigh: {:.3f} [℃], mdot: {:.3f}[kg/s]'.format(Tlow-273.15,Thigh-273.15, self.OutEvap_REF.m))
        print('Tcomp_in: {:.3f} [℃], Tcomp_out: {:.3f} [℃]'.format(self.OutEvap_REF.T-273.15,self.InCond_REF.T-273.15))
        print('Cond_UA: {:.3f} [W/℃], Evap_UA: {:.3f} [W/℃]'.format(outputs.cond_UA, outputs.evap_UA))
        
if __name__ == '__main__':
    DSH_ub = 10.0
    DSC_ub = 10.0
    DSH_lb = 0.1
    DSC_lb = 0.1
    
    fluid_list = ['R1233zd(E)','R1336mzz(Z)','R1224yd(Z)']
    
    num_position = len(fluid_list)*10
    num_time = 20
    
    DSH = np.array([DSH_lb+(DSH_ub-DSH_lb)*np.random.random() for i in range(num_position)])
    DSC = np.array([DSC_lb+(DSC_ub-DSC_lb)*np.random.random() for i in range(num_position)])
    frac_list = np.random.random((num_position,4))
    frac_list = np.array([[f[0]/(f[0]+f[1]+f[2]+f[3]),f[1]/(f[0]+f[1]+f[2]+f[3]),f[2]/(f[0]+f[1]+f[2]+f[3]),f[3]/(f[0]+f[1]+f[2]+f[3])] for f in frac_list])
    ref_list = np.random.randint(0,3,(num_position))
    
    
    xx = np.stack((DSH, DSC, np.squeeze(frac_list[:,1]), np.squeeze(frac_list[:,2]), np.squeeze(frac_list[:,3]),ref_list), axis = 1)
    vv = xx
    
    factor_vec = np.array([[2/(DSH_ub - DSH_lb), 2/(DSC_ub-DSC_lb), 1, 1, 1, 1/(len(fluid_list)-1)]])
    
    w = 1.5*np.repeat(factor_vec, num_position, axis = 0)
    c1 = 3.0*np.repeat(factor_vec, num_position, axis = 0)
    c2 = 4.5*np.repeat(factor_vec, num_position, axis = 0)
    
    opt_result = []
    
    
    for tt in range(num_time):
        print("No.%d flight----------------------------" %tt)
        fit_array = []
        particle_num = -1
        for dsh, dsc, f1, f2, f3, ref_num in xx:
            particle_num = particle_num+1
            
            evapfluid = 'water'
            inevapT = 43.0+273.15
            inevapp = 101300.0
            evapm = 0.0
            InEvap = ProcessFluid(Y={evapfluid:1.0,},m = evapm, T = inevapT, p = inevapp)
            
            outevapT = 38.0+273.15
            outevapp = 101300.0
            OutEvap = ProcessFluid(Y={evapfluid:1.0,},m = evapm, T = outevapT, p = outevapp)
            
            condfluid = 'water'
            incondT = 0.0
            incondp = 0.0
            condm = 0.0
            InCond = ProcessFluid(Y={condfluid:1.0,},m = condm, T = incondT, p = incondp)
            
            outcondT = 108.0 + 273.15
            outcondp = 0.0
            OutCond = ProcessFluid(Y={condfluid:1.0,},m = condm, T = outcondT, p = outcondp)
            
            
            inputs = Settings()
            inputs.second = 'steam'
                
            inputs.cond_dp = 0.01
            inputs.evap_dp = 0.01
            inputs.cond_type = 'phe'
            inputs.evap_type = 'phe'
            inputs.layout = '3eco'
            inputs.cond_T_pp = 2.0
            inputs.evap_T_pp = 2.0
            
            inputs.T_steam = 103.0 + 273.15
            inputs.m_steam = 1.535
            inputs.T_makeup = 95.0 + 273.15
            inputs.m_makeup = inputs.m_steam
            
            
            inputs.eff_LPcompS1 = 0.75
            inputs.eff_LPcompS2 = 0.75
            inputs.eff_HPcompS1 = 0.75
            inputs.eff_HPcompS2 = 0.75
            
            outputs = Outputs()
            
            inputs.DSH = dsh
            inputs.DSC = dsc
            inputs.frac_list = [1-f1-f2-f3, f1, f2, f3]
            inputs.Y = {'REFPROP::'+fluid_list[int(ref_num)]:1.0,}
            try:
                HP1000RT = VCHP(InCond, OutCond, InEvap, OutEvap, inputs)
                (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = HP1000RT()
            except:
                outputs.COP_heating = 0.0
                
            fit_array.append(outputs.COP_heating)
            print("No.%d Particle--> DSH:%.2f[℃], DSC:%.2f[℃], f1:%.2f, f2:%.2f, f3:%.2f, REF:%s, COP:%.2f" %(particle_num, dsh, dsc, f1, f2, f3, fluid_list[int(ref_num)], outputs.COP_heating))
        print("Global best: %.2f (gdx:%d)----------------------------" %(np.max(fit_array),int(np.argmax(fit_array))))
        opt_result.append(np.max(fit_array))
        
        if tt == 0:
            fit_best = np.array(fit_array)
            gdx = int(np.argmax(fit_best))
            g_pos = np.repeat(np.reshape(xx[gdx][:],(1,6)),repeats=num_position, axis=0)
            vv_new = np.multiply(w, vv)+np.multiply(np.multiply(c2,np.random.rand(num_position, 6)), g_pos-xx)
            i_pos = xx
        else:
            idx = np.where(fit_array > fit_best)
            i_pos[idx,:] = xx[idx,:]
            fit_best = np.maximum(fit_array, fit_best)
            gdx = int(np.argmax(fit_best))
            g_pos = np.repeat(np.reshape(xx[gdx][:],(1,6)),repeats=num_position, axis=0)
            vv_new = np.multiply(w, vv)+np.multiply(np.multiply(c1,np.random.rand(num_position, 6)), i_pos-xx)+np.multiply(np.multiply(c2,np.random.rand(num_position, 6)), g_pos-xx)
            
        xx = xx + vv_new
        vv = vv_new
        
        for i in range(num_position):
            xx[i][0] = min(max(DSH_lb, xx[i][0]),DSH_ub)
            xx[i][1] = min(max(DSC_lb, xx[i][1]),DSC_ub)
            xx[i][2] = min(max(0, xx[i][2]),1)
            xx[i][3] = min(max(0, xx[i][3]),1)
            xx[i][4] = min(max(0, xx[i][4]),1)
            xx[i][5] = min(max(0, int(xx[i][5])),2)
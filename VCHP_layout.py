from copy import deepcopy
import HX_module as HX
import COMPAND_module as CP
from CoolProp.CoolProp import PropsSI, set_config_string, ALTERNATIVE_REFPROP_PATH
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
        
        if self.inputs.layout == '3eco':
            (self.InCond, self.OutCond, self.InEvap, self.OutEvap, self.InCond_REF, self.OutCond_REF, self.InEvap_REF, self.OutEvap_REF, outputs) = self.Three_Economizer_Solver(self.InCond, self.OutCond, self.InEvap, self.OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, self.inputs, outputs, no_input, cond_ph, evap_ph)
        
        #self.Post_Processing(self.inputs, outputs)
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
                
            if inputs.layout == '3eco':
                eco1_p = InEvap_REF.p*(OutCond_REF.p/InEvap_REF.p)**self.eco1_frac
                eco1_h_vap = PropsSI("H","P",eco1_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco1_h_liq = PropsSI("H","P",eco1_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                eco2_p = eco1_p*(OutCond_REF.p/InEvap_REF.p)**self.eco2_frac
                eco2_h_vap = PropsSI("H","P",eco2_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco2_h_liq = PropsSI("H","P",eco2_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                eco3_p = eco2_p*(OutCond_REF.p/InEvap_REF.p)**self.eco3_frac    
                eco3_h_vap = PropsSI("H","P",eco3_p,"Q",1.0, OutCond_REF.fluidmixture)
                eco3_h_liq = PropsSI("H","P",eco3_p,"Q",0.0, OutCond_REF.fluidmixture)
                
                outputs.eco3_x = (OutCond_REF.h - eco3_h_liq)/(eco3_h_vap - eco3_h_liq)
                outputs.eco2_x = (eco3_h_liq-eco2_h_liq)/(eco2_h_vap - eco2_h_liq)
                outputs.eco1_x = (eco2_h_liq-eco1_h_liq)/(eco1_h_vap - eco1_h_liq)
                        
                outputs.Out_LPcompS1 = deepcopy(OutEvap_REF)
                outputs.Out_LPcompS1.p = eco1_p
                
                LPcompS1 = CP.Compander_module(OutEvap_REF, outputs.Out_LPcompS1)
                (inputs.DSH_1, cond_a) = LPcompS1.COMP(eff_isen = inputs.eff_LPcompS1, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                outputs.Out_LPcompS1 = LPcompS1.primary_out
                
                outputs.In_LPcompS2 = deepcopy(outputs.Out_LPcompS1)
                outputs.In_LPcompS2.h = outputs.Out_LPcompS1.h*(outputs.eco1_x*(1-inputs.inj_ratio_g)+(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))+eco1_h_vap*outputs.eco1_x*inputs.inj_ratio_g+eco1_h_liq*(1-outputs.eco1_x)*inputs.inj_ratio_l
                outputs.In_LPcompS2.T = PropsSI("T","H",outputs.In_LPcompS2.h,"P",outputs.In_LPcompS2.p,outputs.In_LPcompS2.fluidmixture)
                outputs.In_LPcompS2.s = PropsSI("S","T",outputs.In_LPcompS2.T,"P",outputs.In_LPcompS2.p,outputs.In_LPcompS2.fluidmixture)
                
                outputs.Out_LPcompS2 = deepcopy(outputs.In_LPcompS2)
                outputs.Out_LPcompS2.p = eco2_p
                LPcompS2 = CP.Compander_module(outputs.In_LPcompS2, outputs.Out_LPcompS2)
                (inputs.DSH_2, cond_a) = LPcompS2.COMP(eff_isen = inputs.eff_LPcompS2, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                outputs.Out_LPcompS2 = LPcompS2.primary_out
                                    
                outputs.In_HPcompS1 = deepcopy(outputs.Out_LPcompS2)
                outputs.In_HPcompS1.h = outputs.Out_LPcompS2.h*(outputs.eco2_x*(1-inputs.inj_ratio_g)+(1-outputs.eco2_x)*(1-inputs.inj_ratio_l)) + eco2_h_vap*outputs.eco2_x*inputs.inj_ratio_g+eco2_h_liq*(1-outputs.eco2_x)*inputs.inj_ratio_l
                outputs.In_HPcompS1.T = PropsSI("T","H",outputs.In_HPcompS1.h,"P",outputs.In_HPcompS1.p,outputs.In_HPcompS1.fluidmixture)
                outputs.In_HPcompS1.s = PropsSI("S","T",outputs.In_HPcompS1.T,"P",outputs.In_HPcompS1.p,outputs.In_HPcompS1.fluidmixture)
                
                outputs.Out_HPcompS1 = deepcopy(outputs.In_HPcompS1)
                outputs.Out_HPcompS1.p = eco3_p
                HPcompS1 = CP.Compander_module(outputs.In_HPcompS1, outputs.Out_HPcompS1)
                (inputs.DSH_3, cond_a) = HPcompS1.COMP(eff_isen = inputs.eff_HPcompS1, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                outputs.Out_HPcompS1 = HPcompS1.primary_out
                
                outputs.In_HPcompS2 = deepcopy(outputs.Out_HPcompS1)
                outputs.In_HPcompS2.h = outputs.Out_HPcompS1.h*(outputs.eco3_x*(1-inputs.inj_ratio_g)+(1-outputs.eco3_x)*(1-inputs.inj_ratio_l)) + eco3_h_vap*outputs.eco3_x*inputs.inj_ratio_g+eco3_h_liq*(1-outputs.eco3_x)*inputs.inj_ratio_l
                outputs.In_HPcompS2.T = PropsSI("T","H",outputs.In_HPcompS2.h,"P",outputs.In_HPcompS2.p, outputs.In_HPcompS2.fluidmixture)
                outputs.In_HPcompS2.s = PropsSI("S","T",outputs.In_HPcompS2.T,"P",outputs.In_HPcompS2.p, outputs.In_HPcompS2.fluidmixture)
                
                HPcompS2 = CP.Compander_module(outputs.In_HPcompS2, InCond_REF)
                (inputs.DSH_4, cond_a) = HPcompS2.COMP(eff_isen = inputs.eff_HPcompS2, eff_mech = inputs.mech_eff, DSH = inputs.DSH)
                InCond_REF = HPcompS2.primary_out
                                
                InEvap_REF.h = (eco1_h_vap*outputs.eco1_x*(1-inputs.inj_ratio_g)+eco1_h_liq*(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))/(outputs.eco1_x*(1-inputs.inj_ratio_g)+(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))
                InEvap_REF.T = PropsSI("T","H",InEvap_REF.h,"P",InEvap_REF.p,InEvap_REF.fluidmixture)
                
                
            if (no_input == 'InCondT') or (no_input == 'OutCondT') or (no_input == 'Condm'):
                InEvap_REF.m = InEvap.q/(InEvap_REF.h - OutEvap_REF.h)
                OutEvap_REF.m = InEvap_REF.m
                
                if inputs.layout == '3eco':
                    outputs.Out_LPcompS1.m = OutEvap_REF.m
                    outputs.In_LPcompS2.m = OutEvap_REF.m/(outputs.eco1_x*(1-inputs.inj_ratio_g)+(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))
                    outputs.Out_LPcompS2.m = outputs.In_LPcompS1.m
                    outputs.In_HPcompS1.m = outputs.Out_LPcompS2.m/(outputs.eco2_x*(1-inputs.inj_ratio_g)+(1-outputs.eco2_x)*(1-inputs.inj_ratio_l))
                    outputs.Out_HPcompS1.m = outputs.In_HPcompS1.m
                    outputs.In_HPcompS2.m = outputs.Out_HPcompS1.m/(outputs.eco3_x*(1-inputs.inj_ratio_g)+(1-outputs.eco3_x)*(1-inputs.inj_ratio_l))
                    InCond_REF.m = outputs.In_HPcompS2.m
                    OutCond_REF.m = InCond_REF.m
                    outputs.W_LPcompS1 = LPcompS1.Pspecific*outputs.Out_LPcompS1.m
                    outputs.W_LPcompS2 = LPcompS2.Pspecific*outputs.Out_LPcompS2.m
                    outputs.W_HPcompS1 = HPcompS1.Pspecific*outputs.Out_HPcompS1.m
                    outputs.W_HPcompS2 = HPcompS2.Pspecific*InCond_REF.m
                    
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
                
                if inputs.layout == '3eco':
                    outputs.In_HPcompS2.m = OutCond_REF.m
                    outputs.Out_HPcompS1.m = outputs.In_HPcompS2.m*(outputs.eco3_x*(1-inputs.inj_ratio_g)+(1-outputs.eco3_x)*(1-inputs.inj_ratio_l))
                    outputs.In_HPcompS1.m = outputs.Out_HPcompS1.m
                    outputs.Out_LPcompS2.m = outputs.In_HPcompS1.m*(outputs.eco2_x*(1-inputs.inj_ratio_g)+(1-outputs.eco2_x)*(1-inputs.inj_ratio_l))
                    outputs.In_LPcompS2.m = outputs.Out_LPcompS2.m
                    outputs.Out_LPcompS1.m = outputs.In_LPcompS2.m*(outputs.eco1_x*(1-inputs.inj_ratio_g)+(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))
                    OutEvap_REF.m = outputs.Out_LPcompS1.m
                    InEvap_REF.m = OutEvap_REF.m
                    outputs.W_LPcompS1 = LPcompS1.Pspecific*outputs.Out_LPcompS1.m
                    outputs.W_LPcompS2 = LPcompS2.Pspecific*outputs.Out_LPcompS2.m
                    outputs.W_HPcompS1 = HPcompS1.Pspecific*outputs.Out_HPcompS1.m
                    outputs.W_HPcompS2 = HPcompS2.Pspecific*InCond_REF.m
                    
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
            
            if inputs.layout == '3eco':
                outputs.Wcomp = (outputs.W_LPcompS1+outputs.W_LPcompS2+outputs.W_HPcompS1+outputs.W_HPcompS2)
                outputs.COP_heating = abs(OutCond.q)/outputs.Wcomp
                
        
        outputs.DSH = inputs.DSH
        outputs.evap_UA = evap.UA
        
        return (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    
    def Three_Economizer_Solver(self,InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph, evap_ph):
        self.eco1_frac = inputs.frac_list[0]
        self.eco2_frac = inputs.frac_list[1]
        self.eco3_frac = inputs.frac_list[2]
        
        (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = self.Cycle_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs, no_input, cond_ph, evap_ph)
        
        return(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs)
    
    def Plot_diagram(self, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs):
        (p_array, h_array, T_array, s_array, p_points, h_points, s_points, T_points) = self.Dome_Draw(InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs)
        fig_ph, ax_ph = PLT.subplots()
        p_scale = 1.0e5
        p_unit = '[bar]'
        if inputs.layout == '3eco':
            from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
            ax_ph.set_yscale('log')
            
        ax_ph.plot([i/1.0e3 for i in h_array], [i/p_scale for i in p_array],'k--')
        ax_ph.set_xlabel('Enthalpy [kJ/kg]',fontsize = 15)
        ax_ph.set_ylabel('Pressure '+p_unit,fontsize = 15)
        ax_ph.set_title('Pressure-Enthalpy Diagram\nRefrigerant:{}'.format(list(inputs.Y.keys())[0]),fontsize = 18)
        ax_ph.yaxis.set_major_formatter(ScalarFormatter())
        ax_ph.yaxis.set_minor_formatter(FormatStrFormatter('%.1f'))
        ax_ph.tick_params(axis = 'x', labelsize = 10)
        ax_ph.tick_params(axis = 'y', labelsize = 10, which = 'major')
        ax_ph.tick_params(axis = 'y', labelsize = 9, which = 'minor')
        
        fig_ts, ax_ts = PLT.subplots()
        ax_ts.plot([i/1.0e3 for i in s_array], [i-273.15 for i in T_array],'k--')
        ax_ts.set_xlabel('Entropy [kJ/kg-K]',fontsize = 15)
        ax_ts.set_ylabel('Temperature [℃]',fontsize = 15)
        ax_ts.set_title('Temperature-Entropy Diagram\nRefrigerant:{}'.format(list(inputs.Y.keys())[0]),fontsize = 18)
        ax_ts.tick_params(axis = 'x', labelsize = 10)
        ax_ts.tick_params(axis = 'y', labelsize = 10)
        
        
        if inputs.layout == '3eco':
            ax_ts.plot([i/1.0e3 for i in s_points[1]], [i-273.15 for i in T_points[1]], 'ko-')
            ax_ph.plot([i/1.0e3 for i in h_points[1]], [i/p_scale for i in p_points[1]], 'ko-')
            ax_ts.plot([i/1.0e3 for i in s_points[2]], [i-273.15 for i in T_points[2]], 'ko-')
            ax_ph.plot([i/1.0e3 for i in h_points[2]], [i/p_scale for i in p_points[2]], 'ko-')
            ax_ts.plot([i/1.0e3 for i in s_points[3]], [i-273.15 for i in T_points[3]], 'ko-')
            ax_ph.plot([i/1.0e3 for i in h_points[3]], [i/p_scale for i in p_points[3]], 'ko-')
            ax_ts.plot([i/1.0e3 for i in s_points[4]], [i-273.15 for i in T_points[4]], 'ko-')
            ax_ph.plot([i/1.0e3 for i in h_points[4]], [i/p_scale for i in p_points[4]], 'ko-')
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
        
        
        if inputs.layout == '3eco':
            s_points = dict()
            T_points = dict()
            h_points = dict()
            p_points = dict()
            
            Out_LPcompS1_svap = PropsSI("S","P",outputs.Out_LPcompS1.p, "Q", 1.0, outputs.Out_LPcompS1.fluidmixture)
            Out_LPcompS1_hvap = PropsSI("H","P",outputs.Out_LPcompS1.p, "Q", 1.0, outputs.Out_LPcompS1.fluidmixture)
            Out_LPcompS1_Tvap = PropsSI("T","P",outputs.Out_LPcompS1.p, "Q", 1.0, outputs.Out_LPcompS1.fluidmixture)
            Out_LPcompS1_hliq = PropsSI("H","P",outputs.Out_LPcompS1.p, "Q", 0.0, outputs.Out_LPcompS1.fluidmixture)
            
            eco1_h = (Out_LPcompS1_hvap*(outputs.eco1_x)*(1-inputs.inj_ratio_g)+Out_LPcompS1_hliq*(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))/(outputs.eco1_x*(1-inputs.inj_ratio_g)+(1-outputs.eco1_x)*(1-inputs.inj_ratio_l))
            eco1_s_hp = PropsSI("S","H",eco1_h,"P",outputs.Out_LPcompS1.p, outputs.Out_LPcompS1.fluidmixture)
            
            s_points[1] = [OutEvap_REF_svap, OutEvap_REF.s, outputs.Out_LPcompS1.s, Out_LPcompS1_svap, eco1_s_hp, InEvap_REF.s, OutEvap_REF_svap]
            T_points[1] = [OutEvap_REF_Tvap, OutEvap_REF.T, outputs.Out_LPcompS1.T, Out_LPcompS1_Tvap, Out_LPcompS1_Tvap, InEvap_REF.T, OutEvap_REF_Tvap]
            h_points[1] = [OutEvap_REF.h, outputs.Out_LPcompS1.h, eco1_h, InEvap_REF.h, OutEvap_REF.h]
            p_points[1] = [OutEvap_REF.p, outputs.Out_LPcompS1.p, outputs.Out_LPcompS1.p, InEvap_REF.p, OutEvap_REF.p]
            
            Out_LPcompS2_svap = PropsSI("S","P",outputs.Out_LPcompS2.p, "Q", 1.0, outputs.Out_LPcompS2.fluidmixture)
            Out_LPcompS2_hvap = PropsSI("H","P",outputs.Out_LPcompS2.p, "Q", 1.0, outputs.Out_LPcompS2.fluidmixture)
            Out_LPcompS2_Tvap = PropsSI("T","P",outputs.Out_LPcompS2.p, "Q", 1.0, outputs.Out_LPcompS2.fluidmixture)
            Out_LPcompS2_hliq = PropsSI("H","P",outputs.Out_LPcompS2.p, "Q", 0.0, outputs.Out_LPcompS2.fluidmixture)
            
            eco2_h = (Out_LPcompS2_hvap*(outputs.eco2_x)*(1-inputs.inj_ratio_g)+Out_LPcompS2_hliq*(1-outputs.eco2_x)*(1-inputs.inj_ratio_l))/(outputs.eco2_x*(1-inputs.inj_ratio_g)+(1-outputs.eco2_x)*(1-inputs.inj_ratio_l))
            eco2_s_hp = PropsSI("S","H",eco2_h,"P",outputs.Out_LPcompS2.p,outputs.Out_LPcompS2.fluidmixture)
            eco2_s_lp = PropsSI("S","H",eco2_h,"P",outputs.In_LPcompS2.p,outputs.In_LPcompS2.fluidmixture)
            
            s_points[2] = [Out_LPcompS1_svap, outputs.In_LPcompS2.s, outputs.Out_LPcompS2.s, Out_LPcompS2_svap, eco2_s_hp, eco2_s_lp, Out_LPcompS1_svap]
            T_points[2] = [Out_LPcompS1_Tvap, outputs.In_LPcompS2.T, outputs.Out_LPcompS2.T, Out_LPcompS2_Tvap, Out_LPcompS2_Tvap, Out_LPcompS1_Tvap, Out_LPcompS1_Tvap]
            h_points[2] = [outputs.In_LPcompS2.h, outputs.Out_LPcompS2.h, eco2_h, eco2_h, outputs.In_LPcompS2.h]
            p_points[2] = [outputs.In_LPcompS2.p, outputs.Out_LPcompS2.p, outputs.Out_LPcompS2.p, outputs.In_LPcompS2.p, outputs.In_LPcompS2.p]
            
            Out_HPcompS1_svap = PropsSI("S","P",outputs.Out_HPcompS1.p, "Q", 1.0, outputs.Out_HPcompS1.fluidmixture)
            Out_HPcompS1_hvap = PropsSI("H","P",outputs.Out_HPcompS1.p, "Q", 1.0, outputs.Out_HPcompS1.fluidmixture)
            Out_HPcompS1_Tvap = PropsSI("T","P",outputs.Out_HPcompS1.p, "Q", 1.0, outputs.Out_HPcompS1.fluidmixture)
            Out_HPcompS1_hliq = PropsSI("H","P",outputs.Out_HPcompS1.p, "Q", 0.0, outputs.Out_HPcompS1.fluidmixture)
            
            eco3_h = (Out_HPcompS1_hvap*(outputs.eco3_x)*(1-inputs.inj_ratio_g)+Out_HPcompS1_hliq*(1-outputs.eco3_x)*(1-inputs.inj_ratio_l))/(outputs.eco3_x*(1-inputs.inj_ratio_g)+(1-outputs.eco3_x)*(1-inputs.inj_ratio_l))
            eco3_s_hp = PropsSI("S","H",eco3_h, "P", outputs.Out_HPcompS1.p, outputs.Out_HPcompS1.fluidmixture)
            eco3_s_lp = PropsSI("S","H",eco3_h, "P", outputs.In_HPcompS1.p, outputs.In_HPcompS1.fluidmixture)
            
            s_points[3] = [Out_LPcompS2_svap, outputs.In_HPcompS1.s, outputs.Out_HPcompS1.s, Out_HPcompS1_svap, eco3_s_hp, eco3_s_lp, Out_LPcompS2_svap]
            T_points[3] = [Out_LPcompS2_Tvap, outputs.In_HPcompS1.T, outputs.Out_HPcompS1.T, Out_HPcompS1_Tvap, Out_HPcompS1_Tvap, Out_LPcompS2_Tvap, Out_LPcompS2_Tvap]
            h_points[3] = [outputs.In_HPcompS1.h, outputs.Out_HPcompS1.h, eco3_h, eco3_h, outputs.In_HPcompS1.h]
            p_points[3] = [outputs.In_HPcompS1.p, outputs.Out_HPcompS1.p, outputs.Out_HPcompS1.p, outputs.In_HPcompS1.p, outputs.In_HPcompS1.p]
            
            
            ref_s_lp = PropsSI("S","H",OutCond_REF.h,"P",outputs.In_HPcompS2.p,OutCond_REF.fluidmixture)
            
            s_points[4] = [Out_HPcompS1_svap, outputs.In_HPcompS2.s, InCond_REF.s, InCond_REF_svap, OutCond_REF_sliq, OutCond_REF.s, ref_s_lp, Out_HPcompS1_svap]
            T_points[4] = [Out_HPcompS1_Tvap, outputs.In_HPcompS2.T, InCond_REF.T, InCond_REF_Tvap, OutCond_REF_Tliq, OutCond_REF.T, Out_HPcompS1_Tvap, Out_HPcompS1_Tvap]
            h_points[4] = [outputs.In_HPcompS2.h, InCond_REF.h, OutCond_REF.h, OutCond_REF.h, outputs.In_HPcompS2.h]
            p_points[4] = [outputs.In_HPcompS2.p, InCond_REF.p, OutCond_REF.p, outputs.In_HPcompS2.p, outputs.In_HPcompS2.p]
        return (p_array, h_array, T_array, s_array, p_points, h_points, s_points, T_points)
            
            
    def Post_Processing(self, inputs, outputs):
        import sys
        sys.stdout = open("result.txt","w")
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
        
        if inputs.layout == '3eco':
            print('InEvap_REF (T: %.2f [℃], P: %.2f [bar], H: %.3f [kJ/kg], mdot: %.3f [kg/s])' %(self.InEvap_REF.T-273.15, self.InEvap_REF.p/1.0e5, self.InEvap_REF.h/1.0e3, self.InEvap_REF.m))
            print('OutEvap_REF (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(self.OutEvap_REF.T-273.15, self.OutEvap_REF.p/1.0e5, self.OutEvap_REF.m))
            print('Out_LPcompS1 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.Out_LPcompS1.T-273.15, outputs.Out_LPcompS1.p/1.0e5, outputs.Out_LPcompS1.m))
            print('In_LPcompS2 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.In_LPcompS2.T-273.15, outputs.In_LPcompS2.p/1.0e5, outputs.In_LPcompS2.m))
            print('Out_LPcompS2 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.Out_LPcompS2.T-273.15, outputs.Out_LPcompS2.p/1.0e5, outputs.Out_LPcompS2.m))
            print('In_HPcompS1 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.In_HPcompS1.T-273.15, outputs.In_HPcompS1.p/1.0e5, outputs.In_HPcompS1.m))
            print('Out_HPcompS1 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.Out_HPcompS1.T-273.15, outputs.Out_HPcompS1.p/1.0e5, outputs.Out_HPcompS1.m))
            print('In_HPcompS2 (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(outputs.In_HPcompS2.T-273.15, outputs.In_HPcompS2.p/1.0e5, outputs.In_HPcompS2.m))
            print('InCond_REF (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(self.InCond_REF.T-273.15, self.InCond_REF.p/1.0e5, self.InCond_REF.m))
            print('OutCond_REF (T: %.2f [℃], P: %.2f [bar], mdot: %.3f [kg/s])' %(self.OutCond_REF.T-273.15, self.OutCond_REF.p/1.0e5, self.OutCond_REF.m))
            print('Economizer #3 quality: %.2f' %outputs.eco3_x)
            print('Economizer #2 quality: %.2f' %outputs.eco2_x)
            print('Economizer #1 quality: %.2f' %outputs.eco1_x)
        sys.stdout.close()
if __name__ == '__main__':
    set_config_string(ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
    
    inputs = Settings()
    inputs.second = 'steam'
        
    inputs.cond_dp = 0.005
    inputs.evap_dp = 0.005
    inputs.cond_type = 'phe'
    inputs.evap_type = 'phe'
    inputs.layout = '3eco'
    inputs.cond_T_pp = 5.0
    inputs.evap_T_pp = 5.0
    
    inputs.T_steam = 103.0 + 273.15
    inputs.m_steam = 1.371
    inputs.T_makeup = 30.0 + 273.15
    inputs.m_makeup = inputs.m_steam
    
    inputs.mech_eff = 0.95
    inputs.eff_LPcompS1 = 0.7
    inputs.eff_LPcompS2 = 0.7
    inputs.eff_HPcompS1 = 0.7
    inputs.eff_HPcompS2 = 0.7
    
    inputs.inj_ratio_g = 1.0
    inputs.inj_ratio_l = 0.0
    
    DSH_ub = 20.0
    DSC_ub = 20.0
    DSH_lb = 1.0
    DSC_lb = 1.0
    
    fluid_list = ['R1233zd(E)','R1224yd(Z)']
    
    
    num_position = len(fluid_list)*20
    num_time = 20
    
    DSH = np.array([DSH_lb+(DSH_ub-DSH_lb)*np.random.random() for i in range(num_position)])
    DSC = np.array([DSC_lb+(DSC_ub-DSC_lb)*np.random.random() for i in range(num_position)])
    frac_list = np.random.rand(num_position,4)
    frac_list = [frac_list[i,:]/sum(frac_list[i,:]) for i in range(num_position)]
    frac_list = np.reshape(frac_list,(num_position,4))
    ref_list = np.random.randint(0,len(fluid_list),(num_position))
    
    xx = np.stack((DSH, DSC, np.squeeze(frac_list[:,0]), np.squeeze(frac_list[:,1]), np.squeeze(frac_list[:,2]),ref_list), axis = 1)
    vv = xx
    
    f_coef = 0.1
    
    factor_vec = np.array([[2/(DSH_ub - DSH_lb), 2/(DSC_ub-DSC_lb), f_coef, f_coef, f_coef, 1/(len(fluid_list)-1)]])
    
    w = 1.5*np.repeat(factor_vec, num_position, axis = 0)
    c1 = 3.0*np.repeat(factor_vec, num_position, axis = 0)
    c2 = 4.5*np.repeat(factor_vec, num_position, axis = 0)
    
    opt_result_max = []
    opt_result_avg = []

    
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
            
            outputs = Outputs()
            
            inputs.DSH = dsh
            inputs.DSC = dsc
            inputs.frac_list = [f1, f2, f3]
            inputs.Y = {'REFPROP::'+fluid_list[int(ref_num)]:1.0,}

            HP1000RT = VCHP(InCond, OutCond, InEvap, OutEvap, inputs)
            (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = HP1000RT()

                
            fit_array.append(outputs.COP_heating)
            print("No.%d Particle--> DSH:%.2f[℃], DSC:%.2f[℃], f1:%.2f, f2:%.2f, f3:%.2f, REF:%s, COP:%.2f" %(particle_num, dsh, dsc, f1, f2, f3, fluid_list[int(ref_num)], outputs.COP_heating))
        print("Global best: %.2f (gdx:%d)----------------------------" %(np.max(fit_array),int(np.argmax(fit_array))))
        opt_result_max.append(np.max(fit_array))
        opt_result_avg.append(np.mean(fit_array))
        
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
            xx[i][2] = min(max(0.0, xx[i][2]),1.0)
            xx[i][3] = min(max(0.0, xx[i][3]),1.0)
            xx[i][4] = min(max(0.0, xx[i][4]),1.0)
            xx[i][5] = min(max(0, int(xx[i][5])),len(fluid_list)-1)
        
    fig_max, ax_max = PLT.subplots()
    ax_max.plot(np.arange(1,num_time), opt_result_max[1:], 'kD-')
    ax_max.set_xlabel('Trial Number', fontsize = 15)
    ax_max.set_ylabel('COP', fontsize = 15)
    ax_max.set_title('PSO algorithm results')    
    ax_max.set_xticks(np.arange(0,num_time-1, 5))
    ax_max.tick_params(axis = 'x', labelsize = 13)
    ax_max.tick_params(axis = 'y', labelsize = 13)
    
    fig_avg, ax_avg = PLT.subplots()
    ax_avg.plot(np.arange(1,num_time), opt_result_avg[1:], 'kD-')
    ax_avg.set_xlabel('Trial Number', fontsize = 15)
    ax_avg.set_ylabel('COP', fontsize = 15)
    ax_avg.set_title('PSO algorithm results')
    ax_avg.set_xticks(np.arange(0,num_time-1, 5))
    ax_avg.tick_params(axis = 'x', labelsize = 13)
    ax_avg.tick_params(axis = 'y', labelsize = 13)
    
    fig_max.savefig('.\Figs\Max_diagram.png',dpi=300)
    fig_avg.savefig('.\Figs\Avg_diagram.png',dpi=300)
    # 결과
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
    
    outputs = Outputs()
    
    ggdx = int(np.argmax(fit_array))
    
    inputs.DSH = xx[ggdx][0]
    inputs.DSC = xx[ggdx][1]
    inputs.frac_list = [xx[ggdx][2], xx[ggdx][3], xx[ggdx][4]]
    inputs.Y = {'REFPROP::'+fluid_list[int(xx[ggdx][5])]:1.0,}
    
    HP1000RT = VCHP(InCond, OutCond, InEvap, OutEvap, inputs)
    (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = HP1000RT()
    
    HP1000RT.Post_Processing(inputs, outputs)
    HP1000RT.Plot_diagram(InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs)
    
    
    
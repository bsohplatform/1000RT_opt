from dataclasses import dataclass, field
from re import T
from CoolProp.CoolProp import PropsSI

"""
Common Data Types
"""
@dataclass
class ProcessFluid:  
    def __init__(self,  Y: dict = field(default_factory=dict),  m: float = 0.0,  T: float = 0.0, p: float = 0.0, q: float = 0.0, h: float = 0.0, s: float = 0.0, Cp: float = 0.0):
        self.fluidmixture: str = ''
        if len(list(Y.keys())) == 1:
            self.fluidmixture = list(Y.keys())[0]
        else:
            for fluids, ratio in Y.items():
                if fluids == list(Y.keys())[0]:
                    self.fluidmixture = self.fluidmixture+'SRK::'+fluids+'['+str(ratio)+']'+'&'      
                elif fluids == list(Y.keys())[-1]:
                    self.fluidmixture = self.fluidmixture+fluids+'['+str(ratio)+']'
                else:
                    self.fluidmixture = self.fluidmixture+fluids+'['+str(ratio)+']'+'&'
                        
        self.m = m
        self.T = T
        self.p = p
        self.q = q
        self.h = h
        self.s = s
        self.Cp = Cp
                
        
@dataclass
class Settings:
    # 냉매 입력
    Y = {'R410A':1.0,}
    
    # 공정 정보
    second: str = 'process'
    cycle: str = 'vcc'
    layout: str = 'ihx'
    frac: str = 0.0
    
    DSC = 5.0
    DSH = 10.0
    
    # 응축기 스펙
    cond_type = 'phe'
    cond_T_pp: float = 2.0
    cond_T_lm: float = 10.0
    cond_dp: float = 0.01
    cond_N_element: int = 20
    cond_N_row: int = 5
    cond_UA = 0.0
    
    
    # 증발기 스펙
    evap_type = 'phe'
    evap_T_pp: float = 2.0
    evap_T_lm: float = 10.0
    evap_dp: float = 0.01
    evap_N_element: int = 20
    evap_N_row: int = 5
    evap_UA = 0.0
    
    # 터보기기 스펙
    comp_eff: float = 0.7
    expand_eff: float = 0.0
    mech_eff: float = 1.0
    
    # 증기 공정 조건 입력
    T_steam: float = 120.0
    T_steam = T_steam + 273.15
    m_steam: float = 0.1
    T_makeup: float = 30.0
    T_makeup = T_makeup + 273.15
    m_makeup = m_steam
    
    # Three economizer inputs
    frac_list = [0.25, 0.25, 0.25, 0.25]
    eff_LPcompS1 = 0.75
    eff_LPcompS2 = 0.75
    eff_HPcompS1 = 0.75
    eff_HPcompS2 = 0.75
    
    # 수렴오차
    tol: float = 1.0e-3
    
    
@dataclass
class Outputs:
    COP_heating: float = 0.0
    Wcomp: float = 0.0
    Wexpand: float = 0.0
    cond_UA:float = 0.0
    evap_UA:float = 0.0
    DSH: float = 0.0
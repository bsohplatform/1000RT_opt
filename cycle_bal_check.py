from VCHP_layout import VCHP
from HP_dataclass import ProcessFluid, Settings, Outputs

inputs = Settings()
inputs.second = 'steam'
    
inputs.cond_dp = 0.00
inputs.evap_dp = 0.00
inputs.cond_type = 'phe'
inputs.evap_type = 'phe'
inputs.layout = '3eco'
inputs.cond_T_pp = 5.0
inputs.evap_T_pp = 5.0

inputs.T_steam = 103.0 + 273.15
inputs.m_steam = 1.371
inputs.T_makeup = 30.0 + 273.15
inputs.m_makeup = inputs.m_steam

inputs.mech_eff = 1.0
inputs.eff_LPcompS1 = 0.8
inputs.eff_LPcompS2 = 0.8
inputs.eff_HPcompS1 = 0.8
inputs.eff_HPcompS2 = 0.8

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

inputs.DSH = 0.01
inputs.DSC = 1.0
inputs.inj_ratio_g1 = 0.0
inputs.inj_ratio_g2 = 0.0
inputs.inj_ratio_g3 = 0.0

inputs.inj_ratio_l1 = 0.0
inputs.inj_ratio_l2 = 0.0
inputs.inj_ratio_l3 = 0.0

inputs.frac_list = [0.2, 0.2, 0.2, 0.2]
inputs.Y = {'REFPROP::R1233zd(E)':1.0}

HP1000RT = VCHP(InCond, OutCond, InEvap, OutEvap, inputs)
(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = HP1000RT()

HP1000RT.Post_Processing(inputs, outputs)
HP1000RT.Plot_diagram(InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, inputs, outputs)
import numpy as np
import matplotlib as mpl
import scipy #.optimize

from eos import cs2_qcd_fct
from bjorken_ns_solver import init_bjorken_ns_solver, approx_effective_viscosity_ns, better_effective_viscosity_ns

hbarc=0.1973

# Speed of sound (EOS)
def cs2_fct(T_in_fm):
    return cs2_qcd_fct(T_in_fm)
    #return 1./3

#
#eta/s = 0.12,
#zeta/s (=bulk below):
#double B_norm = 0.15;
#double B_width1 = 0.01/hbarc;
#double B_width2 = 0.10/hbarc;
#double Tpeak = 0.160/hbarc;
#double Tdiff = temperature - Tpeak;
#if (Tdiff > 0.) {
#        Tdiff = Tdiff/B_width2;
#} else {
#        Tdiff = Tdiff/B_width1;
#}
#bulk = B_norm*exp(-Tdiff*Tdiff)

def eta_over_s_fct(T_in_fm):
    return 0.12

def zeta_over_s_fct(T_in_fm):

    B_norm = 0.15
    B_width1 = 0.01/hbarc
    B_width2 = 0.10/hbarc
    Tpeak = 0.160/hbarc
    Tdiff = T_in_fm - Tpeak
    if (Tdiff > 0.):
            Tdiff = Tdiff/B_width2
    else:
            Tdiff = Tdiff/B_width1

    zeta_over_s = B_norm*np.exp(-Tdiff*Tdiff)
    
    return zeta_over_s

def combined_visc_fct(T_in_fm):

    return 4./3.*eta_over_s_fct(T_in_fm)+zeta_over_s_fct(T_in_fm)

# Bjorken evolution parameters
T0_list_in_fm=np.array([0.2,0.3,0.4,0.5,0.6])/hbarc
tau0=0.1
tauf=0.8


for T0_in_fm in T0_list_in_fm:

    ns_sol=init_bjorken_ns_solver(T0_in_fm, tau0, cs2_fct, combined_visc_fct_param=None)

    eff_eta_over_s=3./4.*better_effective_viscosity_ns(tau0, tauf, ns_sol, cs2_fct, combined_visc_fct)

    print("Initial temperature=", T0_in_fm*hbarc, ", effective shear viscosity=", eff_eta_over_s)



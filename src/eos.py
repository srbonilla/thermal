import os
import numpy as np
from numba import njit

class EOSvaporcurve:
    """Class for vapor curve from ANEOS."""
    def __init__(self):
        self.A1_T = None
        self.A1_rho_liq = None
        self.A1_rho_vap = None
        self.A1_P_liq = None
        self.A1_P_vap = None
        self.A1_u_liq = None
        self.A1_u_vap = None
        self.A1_s_liq = None
        self.A1_s_vap = None
        self.units = ''

class EOSmeltcurve:
    """Class for melt curve from ANEOS."""
    def __init__(self):
        self.A1_T  = None
        self.A1_rho_liq = None
        self.A1_rho_sol = None
        self.A1_P_liq = None
        self.A1_P_sol = None
        self.A1_u_liq = None
        self.A1_u_sol = None
        self.A1_s_liq = None
        self.A1_s_sol = None
        self.units = ''
        
class EOScriticalpoint:
    """Class for critical point state from the EOS."""
    def __init__(self):
        self.P   = 0
        self.s   = 0  
        self.T   = 0 
        self.rho = 0
        self.u   = 0
        self.units = ''
        
def read_eosvaporcurve(filepath):
    vc = EOSvaporcurve()
    with open(filepath) as f:
        lines = f.readlines()

    vc.units = lines[1][:-1]
    vc.A1_T = np.array(lines[2].split(", ")[:-1], dtype='float')
    vc.A1_rho_liq = np.array(lines[3].split(", ")[:-1], dtype='float')
    vc.A1_rho_vap = np.array(lines[4].split(", ")[:-1], dtype='float')
    vc.A1_P_liq = np.array(lines[5].split(", ")[:-1], dtype='float')
    vc.A1_P_vap = np.array(lines[6].split(", ")[:-1], dtype='float')
    vc.A1_u_liq = np.array(lines[7].split(", ")[:-1], dtype='float')
    vc.A1_u_vap = np.array(lines[8].split(", ")[:-1], dtype='float')
    vc.A1_s_liq = np.array(lines[9].split(", ")[:-1], dtype='float')
    vc.A1_s_vap = np.array(lines[10].split(", ")[:-1], dtype='float')
    
    return vc

def read_eosmeltcurve(filepath):
    mc = EOSmeltcurve()
    with open(filepath) as f:
        lines = f.readlines()

    mc.units = lines[1][:-1]
    mc.A1_T = np.array(lines[2].split(", ")[:-1], dtype='float')
    mc.A1_rho_liq = np.array(lines[3].split(", ")[:-1], dtype='float')
    mc.A1_rho_sol = np.array(lines[4].split(", ")[:-1], dtype='float')
    mc.A1_P_liq = np.array(lines[5].split(", ")[:-1], dtype='float')
    mc.A1_P_sol = np.array(lines[6].split(", ")[:-1], dtype='float')
    mc.A1_u_liq = np.array(lines[7].split(", ")[:-1], dtype='float')
    mc.A1_u_sol = np.array(lines[8].split(", ")[:-1], dtype='float')
    mc.A1_s_liq = np.array(lines[9].split(", ")[:-1], dtype='float')
    mc.A1_s_sol = np.array(lines[10].split(", ")[:-1], dtype='float')
    
    return mc
        
ANEOS_forsterite_cp = EOScriticalpoint()
ANEOS_forsterite_cp.__dict__ = {'P': 0.159304, 's': 0.00637921, 'T': 6043.58, 'rho': 570.535, 'u': 16.9138, 'units': 'T K, rho kg/m3, P GPa, U MJ/kg, S MJ/K/kg'}

ANEOS_Fe85Si15_cp = EOScriticalpoint()
ANEOS_Fe85Si15_cp.__dict__ = {'P': 0.540777, 's': 0.00402756, 'T': 7279.36, 'rho': 857.1319999999999, 'u': 11.0937, 'units': 'T K, rho kg/m3, P GPa, U MJ/kg, S MJ/K/kg'}

ANEOS_forsterite_vc = read_eosvaporcurve("../data/ANEOS_forsterite_vc.txt")
ANEOS_forsterite_mc = read_eosmeltcurve("../data/ANEOS_forsterite_mc.txt")
ANEOS_Fe85Si15_vc = read_eosvaporcurve("../data/ANEOS_Fe85Si15_vc.txt")
ANEOS_Fe85Si15_mc = read_eosmeltcurve("../data/ANEOS_Fe85Si15_mc.txt")
    
    
def plot_ANEOS_rhou(ax, mat, **kwargs):
    
    lw = kwargs.get("lw", 1)
    
    if mat == "ANEOS_forsterite":
        vc = ANEOS_forsterite_vc
        mc = ANEOS_forsterite_mc
    elif mat == "ANEOS_Fe85Si15":
        vc = ANEOS_Fe85Si15_vc
        mc = ANEOS_Fe85Si15_mc
    else:
        raise NotImplementedError("Unknown material")
        
    ax.plot(vc.A1_rho_liq, vc.A1_u_liq, color="black", linewidth=lw)
    ax.plot(vc.A1_rho_vap, vc.A1_u_vap, color="black", linewidth=lw)
    ax.plot(vc.A1_rho_liq[0], vc.A1_u_liq[0],'ko')
    ax.plot(mc.A1_rho_liq, mc.A1_u_liq, color="black", linewidth=lw)
    ax.plot(mc.A1_rho_sol, mc.A1_u_sol, color="black", linewidth=lw)
    
def plot_ANEOS_sP(ax, mat, **kwargs):
    
    lw = kwargs.get("lw", 1)
    
    if mat == "ANEOS_forsterite":
        vc = ANEOS_forsterite_vc
        mc = ANEOS_forsterite_mc
    elif mat == "ANEOS_Fe85Si15":
        vc = ANEOS_Fe85Si15_vc
        mc = ANEOS_Fe85Si15_mc
    else:
        raise NotImplementedError("Unknown material")
        
    ax.plot(vc.A1_s_liq, vc.A1_P_liq, color="black", linewidth=lw)
    ax.plot(vc.A1_s_vap, vc.A1_P_vap, color="black", linewidth=lw)
    ax.plot(vc.A1_s_liq[0], vc.A1_P_liq[0],'ko')
    ax.plot(mc.A1_s_liq, mc.A1_P_liq, color="black", linewidth=lw)
    ax.plot(mc.A1_s_sol, mc.A1_P_sol, color="black", linewidth=lw)
    
@njit()
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

#@njit()
def above_vc(rho, u, mat):
    
    if mat == "ANEOS_forsterite":
        vc = ANEOS_forsterite_vc
        cp = ANEOS_forsterite_cp
    elif mat == "ANEOS_Fe85Si15":
        vc = ANEOS_Fe85Si15_vc
        cp = ANEOS_Fe85Si15_cp
    else:
        raise NotImplementedError("Unknown material")
        
    # check if density is too high
    if rho > np.max(vc.A1_rho_liq):
        return None
    
    if rho > cp.rho:
        idx, val = find_nearest(vc.A1_rho_liq, rho)
        if u > vc.A1_u_liq[idx]:
            return 1
        else:
            return 0
    else:
        idx, val = find_nearest(vc.A1_rho_vap, rho)
        if u > vc.A1_u_vap[idx]:
            return 1
        else:
            return 0
        
#@njit()
def above_mc_liq(rho, u, mat):
    
    if mat == "ANEOS_forsterite":
        mc = ANEOS_forsterite_mc
    elif mat == "ANEOS_Fe85Si15":
        mc = ANEOS_Fe85Si15_mc
    else:
        raise NotImplementedError("Unknown material")
        
    # check if density is too low
    if rho < np.min(mc.A1_rho_liq):
        return None
    
    
    idx, val = find_nearest(mc.A1_rho_liq, rho)
    if u > mc.A1_u_liq[idx]:
        return 1
    else:
        return 0
    
#@njit()
def above_mc_sol(rho, u, mat):
    
    if mat == "ANEOS_forsterite":
        mc = ANEOS_forsterite_mc
    elif mat == "ANEOS_Fe85Si15":
        mc = ANEOS_Fe85Si15_mc
    else:
        raise NotImplementedError("Unknown material")
        
    # check if density is too low
    if rho < np.min(mc.A1_rho_sol):
        return None
    
    
    idx, val = find_nearest(mc.A1_rho_sol, rho)
    if u > mc.A1_u_sol[idx]:
        return 1
    else:
        return 0
    
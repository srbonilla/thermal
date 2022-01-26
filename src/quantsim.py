import numpy as np
import swiftsimio as sio
import src.vars

class SimData:
    def __init__(self, fp, type):
        assert type in ("settling", "impact")
        
        self.type = type
        self.data = sio.load(fp)
        self.fp = fp
        self.A2_pos = self.data.gas.coordinates.value - self.data.metadata.boxsize.value/2 # R_earth units
        self.A1_r = np.sqrt(np.sum((self.A2_pos)**2, axis=1))
        self.A1_rho = self.data.gas.densities.value * src.vars.M_earth / src.vars.R_earth**3 # kg m3
        self.A1_u = self.data.gas.internal_energies.value * src.vars.R_earth**2 / 1e6 # MJ
        self.A1_P = self.data.gas.pressures.value * src.vars.M_earth / src.vars.R_earth # Pa
        
        
    @property
    def r_eq(self):
        if self.type == "settling":
            mask = np.abs(self.A2_pos[:,2]) < np.min(self.data.gas.smoothing_lengths.value)
            r_eq = np.max(self.A1_r[mask])
            return r_eq
        else:
            return None
        
    @property
    def r_po(self):
        if self.type == "settling":
            A1_r_xy = np.sqrt(self.A2_pos[:,0]**2 + self.A2_pos[:,1]**2)
            # polar radius
            mask = np.abs(A1_r_xy) < np.min(self.data.gas.smoothing_lengths.value)
            r_po = np.max(self.A1_r[mask])
            return r_po
        else:
            return None
        
    @property
    def radius(self):
        if self.type == "settling":
            radius = np.mean([self.r_eq, self.r_po])
            return radius
        else:
            return None
        
    @property
    def total_mass(self):
        mass = np.sum(self.data.gas.masses.value) # M_earth units
        return mass
        
        
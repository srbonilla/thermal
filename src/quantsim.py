import numpy as np
import swiftsimio as sio
import src.vars
import woma

class SimData:
    def __init__(self, fp, type):
        assert type in ("settling", "impact")
        
        self.type = type
        self.data = sio.load(fp)
        self.fp = fp
        
        self.A2_pos = self.data.gas.coordinates.value - self.data.metadata.boxsize.value/2 # R_earth units
        self.A2_vel = self.data.gas.velocities.value
        self.A1_rho = self.data.gas.densities.value * src.vars.M_earth / src.vars.R_earth**3 # kg m3
        self.A1_u = self.data.gas.internal_energies.value * src.vars.R_earth**2 / 1e6 # MJ
        self.A1_P = self.data.gas.pressures.value * src.vars.M_earth / src.vars.R_earth # Pa
        self.A1_m = self.data.gas.masses.value
        self.A1_mat_id = self.data.gas.material_ids.value
        self.A1_h = self.data.gas.smoothing_lengths.value
        self.A1_id = self.data.gas.particle_ids.value
        
        self.units = {
            "A2_pos":"R_earth",
            "A2_vel":"R_earth/s",
            "A1_r":"R_earth",
            "A1_rho":"kg/m3",
            "A1_u":"MJ/kg",
            "A1_P":"Pa",
            "A1_m":"M_earth",
            "A1_L":"L_em",
            "A1_h":"R_earth",
            "A1_T":"K",
            "A1_E":"R_earth^2/s2"
        }
        
    @property
    def A1_r(self):
        A1_r = np.sqrt(np.sum((self.A2_pos)**2, axis=1))
        return A1_r
    
    @property
    def A1_v(self):
        A1_v = np.sqrt(np.sum((self.A2_vel)**2, axis=1))
        return A1_v
    
    @property
    def A1_E(self):
        A1_E = self.data.gas.potentials.value + 0.5*np.sum(self.A2_vel**2, axis=1)
        return A1_E
        
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
        mass = np.sum(self.A1_m) # M_earth units
        return mass
    
    @property
    def A2_L(self):
        A2_L = self.A1_m.reshape(-1,1)*np.cross(self.A2_pos, self.A2_vel)
        A2_L *= src.vars.M_earth * src.vars.R_earth**2 / src.vars.L_em
        return A2_L
    
    @property
    def A1_L(self):
        A1_L = np.sqrt(np.sum((self.A2_L)**2, axis=1))
        return A1_L
        
    
    @property
    def angular_momentum(self):
        L = np.sum(self.A1_L)
        return L
    
    @property
    def centre_of_mass(self):
        CoM = np.sum(self.A1_m.reshape(-1,1)*self.A2_pos, axis=0)/np.sum(self.A1_m)
        return CoM
    
    # as method because slow to compute
    def compute_A1_T(self):
        woma.load_eos_tables()
        A1_T = woma.eos.eos.A1_T_u_rho(self.A1_u * 1e6, self.A1_rho, self.A1_mat_id)
        self.A1_T = A1_T
        
    @property
    def N(self):
        N = len(self.A1_m)
        return N
    
    def use_main_planet_reference(self, mat_core="ANEOS_Fe85Si15", R_roche=3):
        
        assert R_roche < 100, f"R_roche must be in R_earth units"
        self.R_roche = R_roche
        
        mat_id = src.vars.Di_mat_id[mat_core]
        planet_centre = np.median(self.A2_pos[self.A1_mat_id == mat_id, :], axis=0)
        
        self.A2_pos -= planet_centre
        
        mask_in_core = np.logical_and.reduce((self.A1_r < R_roche, self.A1_mat_id == mat_id))
        v_in_core = np.median(self.A2_vel[mask_in_core], axis=0)
        
        self.A2_vel -= v_in_core
        
    @property
    def M_in(self):
        
        assert self.R_roche
        
        mask = self.A1_r < self.R_roche
        M_in = np.sum(self.A1_m[mask])
        
        return M_in
        
        
    
        
        
        
        
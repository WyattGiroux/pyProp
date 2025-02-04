from os import getcwd
from os.path import dirname, join, normpath, isabs
import numpy as np
from toml import load
from scipy.interpolate import RegularGridInterpolator

class CompressorMap:
    def __init__(self, name, defaultmap_name=None, usermap_path=None, boundsError=False, method='linear', fill_value=None):
        if defaultmap_name is None and usermap_path is None:
            raise ValueError(f'Either defaultmap_name or usermap_path must be specified')
        
        if defaultmap_name is not None:
            filepath = join(dirname(__file__), '../data/', defaultmap_name)
        else:
            filepath = normpath(join(getcwd(), usermap_path)) if not isabs(usermap_path) else usermap_path
                
        with open(filepath) as mapfile:
            map_raw = load(mapfile)
        
        self.name = name
        
        # ----- LOAD MAP DATA (DO NOT STORE IN CLASS) ----- #
        # Load the design values of the indepenedent map variables    
        self.Nc_des = map_raw['design']['des_Nc']
        self.R_des = map_raw['design']['des_Rline']
        
        # Load independent grid vectors
        Nc_grid = np.array(map_raw['grid']['grid_Nc'])
        R_grid = np.array(map_raw['grid']['grid_Rline'])
        
        # Get design value indices
        Nc_des_ind = np.where(Nc_grid == self.Nc_des)[0][0]
        R_des_ind  = np.where(R_grid  == self.R_des )[0][0]
        
        # Load characteristic values
        mc = np.array(map_raw['characteristic']['massflow'])
        efi = np.array(map_raw['characteristic']['isen_efficiency'])
        pr = np.array(map_raw['characteristic']['pressure_ratio'])
        
        # ----- SAVE CLASS DATA ----- #
        # Get design values of characteristics
        self.mc_des = mc[Nc_des_ind, R_des_ind]
        self.efi_des = efi[Nc_des_ind, R_des_ind]
        self.pr_des = pr[Nc_des_ind, R_des_ind]
        
        # Get map boundary coordinates
        self.R_surge, self.R_windmill = R_grid[0],  R_grid[-1]
        self.N_min,   self.N_max      = Nc_grid[0], Nc_grid[-1]        
        
        # Initialize map scaling as unscaled
        self.s_pr  = 1.0
        self.s_mc  = 1.0
        self.s_Nc  = 1.0
        self.s_efi = 1.0 
        
        # Create Interpolators
        self.mc_itp  = RegularGridInterpolator((Nc_grid, R_grid), mc,  method=method, bounds_error=boundsError, fill_value=None)
        self.pr_itp  = RegularGridInterpolator((Nc_grid, R_grid), pr,  method=method, bounds_error=boundsError, fill_value=None)
        self.efi_itp = RegularGridInterpolator((Nc_grid, R_grid), efi, method=method, bounds_error=boundsError, fill_value=None)
        
        
    def set_map_scaling(self, pr_des, mc_des, Nc_des, eff_des):
        self.s_pr = (pr_des - 1) / (self.pr_des - 1)
        self.s_mc = mc_des / self.mc_des
        self.s_efi = eff_des / self.efi_des
        self.s_Nc = Nc_des / self.Nc_des
        
        
    def surge_margin(self, Nc_descl, mc, pr):
        pr_surge_descl = self.pr_itp((Nc_descl, self.R_surge))
        mc_surge_descl = self.mc_itp((Nc_descl, self.R_surge))
        
        mc_surge = self.s_mc * mc_surge_descl
        pr_surge = 1.0 + self.s_pr * (pr_surge_descl - 1)
        
        return (pr_surge * mc) / (pr * mc_surge)
    
    
    def evaluate_map(self, Nc, R):
        Nc_descl = Nc / self.s_Nc
        
        mc_descl = self.mc_itp((Nc_descl, R))
        pr_descl = self.pr_itp((Nc_descl, R))
        efi_descl = self.efi_itp((Nc_descl, R))
        
        mc = self.s_mc * mc_descl
        pr = 1.0 + self.s_pr * (pr_descl - 1)
        efi = self.s_efi * efi_descl
        
        sm = self.surge_margin(Nc_descl, mc, pr)
        
        return mc, pr, efi, sm

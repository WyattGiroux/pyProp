from os.path import dirname, join
import numpy as np
from toml import load

class CompressorMap:
    def __init__(self, mapname, relpath='../data/'):
        with open(join(dirname(__file__), relpath, mapname)) as mapfile:
            map_raw = load(mapfile)
        
        # Load the design values of the indepenedent map variables    
        self.Nc_des = map_raw['design']['des_Nc']
        self.R_des = map_raw['design']['des_Rline']
        
        # Load independent grid vectors
        self.Nc_grid = np.array(map_raw['grid']['grid_Nc'])
        self.R_grid = np.array(map_raw['grid']['grid_Rline'])
        
        # Get design value indices
        self.Nc_des_ind = np.where(self.Nc_grid == self.Nc_des)[0]
        self.R_des_ind  = np.where(self.R_grid  == self.R_des )[0]
        
        # Load characteristic values
        self.mc = np.array(map_raw['characteristic']['massflow'])
        self.efi = np.array(map_raw['characteristic']['isen_efficiency'])
        self.pr = np.array(map_raw['characteristic']['pressure_ratio'])
        
        # Get design values of characteristics
        self.mc_des = self.mc[self.Nc_des_ind, self.R_des_ind]
        self.efi_des = self.efi[self.Nc_des_ind, self.R_des_ind]
        self.pr_des = self.pr[self.Nc_des_ind, self.R_des_ind]
        
        # Initialize map scaling as unscaled
        self.s_pr  = 1.0
        self.s_mc  = 1.0
        self.s_Nc  = 1.0
        self.s_efi = 1.0 
        
        
    def set_map_scaling(self, pr_des, mc_des, Nc_des, eff_des):
        self.s_pr = (pr_des - 1) / (self.pr_des - 1)
        self.s_mc = mc_des / self.mc_des
        self.s_efi = eff_des / self.efi_des
        self.s_Nc = Nc_des / self.Nc_des
        
        
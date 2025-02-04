from ..structs.element import Element
from ..structs.option import Option
from ..structs.ports.flowstation import FlowStation
from ..structs.ports.shaftport import ShaftPort

from ..structs.compressor_map import CompressorMap

from ..data.constants import TSL, PSL

import numpy as np

class Compressor(Element):
    def __init__(self, name, **kwargs):
        super().__init__(name)
        
        ###############################
        #    COMPRESSOR ATTRIBUTES    #
        ###############################
        # Operation Speed Variables
        self.Nc = 0
        self.NcDes = 0
        self.NcPct = 0
        self.NcqNcDes = 0
        self.Nmech = 0
        
        # Pressure/Temperature Ratio Variables
        self.PR = 1
        self.PRdes = kwargs.get('PRdes', 1)
        self.TR = 1
        self.TRdes = 1
        
        # Map Data
        self.__loadMap(**kwargs)
        self.Rdes = self.map.R_des
        self.R = self.Rdes
        
        # Heat Transfer
        # TODO: Not implemented
        
        # Stall Margins
        self.SMN = 0
        self.SMM = 0
        
        # Physical Flow
        self.mc = 0
        self.mcDes = 0
        self.mcqmcDes = 0
        
        # Efficiency
        self.eff = 0
        self.effDes = kwargs.get('effDes', 1)
        self.effPoly = 0
        
        # Bleed Flows
        # TODO: Not implemented
        
        # Shaft and Power
        self.pwr = 0
        self.pwrBld = 0
        self.trq = 0
        
        # Inlet and Outlet FlowStations
        self.ports['Fl_I'] = FlowStation('Fl_I')
        self.ports['Fl_O'] = FlowStation('Fl_O')
        
        # Shaft Port
        self.ports['Sh_O'] = ShaftPort('Sh_O')
        
        # On/Off Design Option
        self.options['switchDes'] = Option(
            description='Specifies on or off design operation',
            trigger=True,
            allowedValues=('design', 'offdesign'),
            default='design'
        )
    
    
    def __loadMap(self, **kwargs):
        map_name = kwargs.get('map_name', self.name)
        map_default_name = kwargs.get('map_default_name', None)
        map_user_path = kwargs.get('map_user_path', None) # Assume map data is in the pyProp/data directory
        map_boundserror = kwargs.get('map_boundserror', False) # Allow map extrapolation by default
        map_method = kwargs.get('map_method', 'linear') # Default to linear extrapolation
        map_fill_value = kwargs.get('map_fill', None)
        
        if map_method not in ['linear', 'slinear', 'cubic', 'quintic']:
            raise ValueError(f'{map_method} is not a valid method for map interp/extrapolation. \
                Use:\n> linear\n> slinear\n> cubic\n> quintic')
        
        self.map = CompressorMap(map_name, defaultmap_name=map_default_name, 
                                 usermap_path=map_user_path, boundsError=map_boundserror, 
                                 method=map_method, fill_value=map_fill_value)
    
    
    def runelement(self):
        self.preexecute()
        
        # Temperature and pressure relative to standard day sea level
        thetaIn = self.Fl_I.Tt / TSL
        deltaIn = self.Fl_I.Pt / PSL
        
        # Get corrected mass flow at inlet
        self.mc = self.Fl_I.mdot * np.sqrt(thetaIn) / deltaIn
        
        # Get shaft speed + calculate corrected speed
        self.Nmech = self.Sh_O.Nmech
        self.Nc = self.Nmech / np.sqrt(thetaIn)
        
        if self.options['switchDes'].state == 'design':
            # Get on-design properties set by user
            self.PR = self.PRdes
            self.eff = self.effDes
            
            # Set design mass flow
            self.mcDes = self.mc
            self.mcqmcDes = self.mc / self.mcDes
            
            # Set design corrected speed
            self.NcDes = self.Nc
            
            self.map.set_map_scaling(self.PRdes, self.mcDes, self.NcDes, self.effDes)
        else:
            self.mc, self.PR, self.eff, self.SMN = self.map.evaluate_map(self.Nc, self.R)
        
        # Calculate fractional corr. speed and percent speed
        self.NcqNcDes = self.Nc / self.NcDes
        self.NcPct = self.NcqNcDes * 100
        
        # Compute the ideal (isentropic) exit properties
        idealFlO = self.Fl_I.copy()
        
        sIdeal = self.Fl_I.st
        Ptout = self.Fl_I.Pt * self.PR
        
        idealFlO.setTotal_sP(sIdeal, Ptout)
        
        # Use isentropic efficiency to calculate actual exit properties
        htout = (idealFlO.ht - self.Fl_I.ht) / self.eff + self.Fl_I.ht
        # print(self.Fl_I.ht)
        self.Fl_O.setTotal_hP(htout, Ptout)
        
        # Calculate temperature ratio across compressor
        self.TR = self.Fl_O.Tt / self.Fl_I.Tt
        if self.options['switchDes'].state == 'design':
            self.TRdes = self.TR
            
        # Calculate polytropic efficieny from the Gibbs eq. definition
        ds = self.Fl_O.st - self.Fl_I.st

        self.effPoly = (self.Fl_I.Rt * np.log(self.PR)) / (self.Fl_I.Rt * np.log(self.PR) + ds)
        
        # Calculate required power (define as positive)
        self.pwr = self.Fl_I.mdot * (self.Fl_O.ht - self.Fl_I.ht)
        
        self.calculate()
        
        self.postexecute()
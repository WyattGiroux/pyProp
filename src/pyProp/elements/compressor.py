from ..structs.element import Element
from ..structs.option import Option
from ..structs.ports.flowstation import FlowStation
from ..structs.ports.shaftport import ShaftPort

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
        # TODO: Load map logic
        self.map = kwargs.get('mapName', None)
        self.s_pr  = 1
        self.s_mc  = 1
        self.s_eff = 1
        self.s_Nc  = 1
        
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
            allowedValues=('design', 'offdesign')
        )
    
    
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
        
        # Calculate fractional corr. speed and percent speed
        self.NcqNcDes = self.Nc / self.NcDes
        self.NcPct = self.NcqNcDes * 100
        
        tempFlO = self.Fl_I.copy()       
        
        self.calculate()
        
        self.postexecute()
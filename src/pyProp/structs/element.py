class Element:
    def __init__(self, name):
        self.name = name
        self.ports = {}
        self.options = {}
        
    def preexecute(self):
        pass
    
    
    def calculate(self):
        pass
    
    
    def postexecute(self):
        pass
    
    
    def runelement(self):
        self.preexecute()
        self.calculate()
        self.postexecute()
        
        
    def __getattr__(self, name):
        if name in self.ports.keys():
            return self.ports[name]
        else:
            raise ValueError(f"Unrecognized variable {name}")
        
        
    def setOption(self, option, value):
        if option not in self.ports.keys():
            raise ValueError(f"Unrecognized option {option}")
        if value not in self.options[option].allowedValues:
            raise ValueError(f"Invalid value ({value}) for option variable {option}")
        self.options[option].curState = value
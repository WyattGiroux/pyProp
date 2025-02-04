from .port import Port

class ShaftPort(Port):
    '''Port used for rigid shaft connections
    
    '''
    def __init__(self, name):
        super().__init__(name)
        self.Nmech = 0
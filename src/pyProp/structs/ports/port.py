class Port:
    '''Port parent class
    
    Attributes:
        name (str): Name of the port.
    
    '''
    __slots__ = ('name',)
    
    def __init__(self, name):
        self.name = name
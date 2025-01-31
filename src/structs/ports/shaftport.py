from .port import Port

class ShaftPort(Port):
    def __init__(self, name):
        super().__init__(name)
        self.Nmech = 0
class Option:
    def __init__(self, description, trigger, allowedValues, default=0):
        self.description = description
        self.trigger = trigger
        self.allowedValues = allowedValues
        
        self.state = default
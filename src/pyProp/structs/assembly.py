import networkx as nx

from .element import Element
from .option import Option

class Assembly:
    def __init__(self, name):
        self.name = name
        self.elements = {}
        self.options = {}
        self.connectGraph = nx.Graph()
        
    
    def add_element(self, element:Element):
        self.elements[element.name] = element
        self.connectGraph.add_node(element.name)
        
    
    def remove_element(self, element:str):
        if element not in self.elements.keys():
            raise ValueError(f"{element} does not exist in the Assembly")
        self.connectGraph.remove_node(element)
        del self.elements[element]
        
        
    def add_linkage(self, eName1:str, port1:str, eName2:str, port2:str, linkName):
        # Check if the specified elements have been added to the assembly
        if eName1 not in self.elements.keys():
            raise ValueError(f"The {eName1} element has not been added to the Assembly")
        if eName2 not in self.elements.keys():
            raise ValueError(f"The {eName2} element has not been added to the Assembly")
        
        # Check if the specified ports exist in the associated elements
        if port1 not in self.elements[eName1].ports:
            raise ValueError(f"The {port1} port does not exist in {eName1}")
        if port2 not in self.elements[eName2].ports:
            raise ValueError(f"The {port2} port does not exist in {eName2}")
        
        # Check if the linkage name has already been used
        duplicateNames = [(u, v) for u, v, d in self.connectGraph.edges(data=True) if d.get('name') == linkName]
        if duplicateNames == []:        
            self.connectGraph.add_edge(eName1, eName2, name=linkName, p1=port1, p2=port2)
        else:
            raise ValueError(f"The linkage {linkName} already exists")
        
        
    def remove_linkage(self, linkName):
        edges_to_remove = [(u, v) for u, v, d in self.connectGraph.edges(data=True) if d.get('name') == linkName]
        self.connectGraph.remove_edge(*edges_to_remove)
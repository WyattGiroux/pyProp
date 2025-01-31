'''``pyProp`` uses assembly objects to store elements of a propulsion system, their linkages,
and any option variables needed.
'''
import networkx as nx

from .element import Element
from .option import Option

class Assembly:
    ''' ``pyProp`` Assembly Object

    Attributes:
        name (string): Assembly name
        elements (dict): Dictionary containing all elements in an assembly.
            Keys are the element names, values are child classes of Element() (see
            :ref:`elements`).
        options (dict): Dictionary containing all option variables for an assembly.
            Keys are option names, values are instances of Option().
        connectGraph (networkx.Graph): Graph of elements and their linkages.
    '''
    def __init__(self, name):
        ''' Assembly.__init__()
        
        Args:
            name (str): Name of the assembly.
        '''
        self.name = name
        self.elements = {}
        self.options = {}
        self.connectGraph = nx.Graph()
        
    
    def add_element(self, element:Element):
        '''Adds an element to the assembly.
        
        Args:
            element (Element): Element to be added.
        '''
        self.elements[element.name] = element
        self.connectGraph.add_node(element.name)
        
    
    def remove_element(self, element:str):
        '''Removes an element from an assembly.
        
        Args:
            element (str): Name of the element to be removed.
        '''
        if element not in self.elements.keys():
            raise ValueError(f"{element} does not exist in the Assembly")
        self.connectGraph.remove_node(element)
        del self.elements[element]
        
        
    def add_linkage(self, eName1:str, port1:str, eName2:str, port2:str, linkName):
        '''Connects two elements in an assembly.
        
        Args:
            eName1 (str): Name of the first element.
            port1 (str): Port in the first element to be connected.
            eName2 (str): Name of the second element.
            port2 (str): Port in the second element to be connected.
            linkName (str): Name of the linkage.
        '''
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
        '''Removes a linkage between two elements.
        
        Args:
            linkName (str): Name of the linkage to be removed.
        '''
        edges_to_remove = [(u, v) for u, v, d in self.connectGraph.edges(data=True) if d.get('name') == linkName]
        self.connectGraph.remove_edge(*edges_to_remove)
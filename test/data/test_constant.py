############################################################
#                                                          #
#               Tests For Useful Constants                 #
#               Author: Wyatt Giroux                       #
#               Date: 1/10/25                              #
#                                                          #
############################################################
from pytest import approx, raises
from src.data.constants import *

class TestConstants:
    def test_Rmol(self):
        assert Rmol == approx(8.3144)
    
    def test_gee(self):
        assert gee == approx(9.80665)
        
    def test_βT(self):
        assert βT == approx(-0.0065)
        
    def test_TSL(self):
        assert TSL == approx(288.15)
        
    def test_PSL(self):
        assert PSL == approx(101325)
        
    def test_ρSL(self):
        assert ρSL == approx(1.225)
    
    def test_aSL(self):
        assert aSL == approx(340.294)
        
    def test_Hptrop(self):
        assert Hptrop == approx(11000)
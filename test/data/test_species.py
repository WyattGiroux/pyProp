############################################################
#                                                          #
#               Tests For NASA-9 Species                   #
#               Author: Wyatt Giroux                       #
#               Date: 1/10/25                              #
#                                                          #
############################################################
from pytest import approx, raises
from src.data.species import *

class TestSpecies:
    # Only test a selection of properties to check that .inp file
    # has probably been read correctly
    def test_air_mw(self):
        assert Air.molecweight == approx(28.9651159)
        
    def test_Ar_range0_coeffs(self):
        assert Ar.ranges[(200.0, 1000.0)]['coeffs'] == approx([0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967491])
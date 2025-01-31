from pytest import approx, raises
from src.utils.compressible import *

class TestConstants:
    def test_cp(self):
        assert cp(1.4, 287) == approx(1004.5)
        
    def test_cv(self):
        assert cv(1.4, 287) == approx(717.5)
        
        
class TestIsentropic:
    # Pressure <- Mach Relation
    def test_pqPt(self):
        assert pqPt(1) == approx(0.528281787717)
        
    def test_pqPt_err(self):
        with raises(ValueError):
            pqPt(-1)
            
    # Temperature <- Mach Relation
    def test_TqTt(self):
        assert TqTt(1) == approx(0.83333333333)
        
    def test_TqTt_err(self):
        with raises(ValueError):
            TqTt(-1)
            
    # Mach <- Pressure Relation
    def test_MfromP(self):
        assert MfromP(1.9) == approx(1.00319241054)
        
    def test_MfromP_err(self):
        with raises(ValueError):
            MfromP(0.9)
            
    # Mach <- Temperature Relation
    def test_MfromT(self):
        assert MfromT(1.2) == approx(1.)
        
    def test_MfromT_err(self):
        with raises(ValueError):
            MfromT(0.9)
            

class TestCMF:
    def test_cmf(self):
        assert cmf(100, 101325, 288.15, 1) == approx(0.0167530042838)
        
    def test_Dm(self):
        assert Dm(1) == approx(0.0404184198941)
        
    def test_MfromD_subsonic(self):
        assert MfromD(Dm(0.5), supersonic=False) == approx(0.5)
    
    def test_MfromD_supersonic(self):
        assert MfromD(Dm(1.5), supersonic=True) == approx(1.5)
        
    def test_MfromD_sonicLeft(self):
        assert MfromD(Dm(1.0), supersonic=False) == approx(1, 1e-5)
        
    def test_MfromD_sonicRight(self):
        assert MfromD(Dm(1.0), supersonic=True) == approx(1, 1e-5)
        
    def test_MfromD_err(self):
        with raises(ValueError):
            MfromD(Dm(1.0) * 1.05, supersonic=False)
            

class TestBasics:
    def test_geth(self):
        assert geth(288) == approx(289296)
        
    def test_ds(self):
        assert ds(1000000, 1000) == approx(0.590153499252)
        
    def test_ds_err(self):
        with raises(ValueError):
            ds(-1, 1)
            ds(1, -1)
    
    def test_impulse(self):
        assert impulse(101325, 1, 100, 135.343) == approx(114859.3)
    
    def test_htrArea(self):
        assert htrArea(1, .5) == approx(0.589048622548)
        
    def test_htrAreaMean(self):
        assert htrAreaMean(0.375, 0.5) == approx(0.589048622548)
        
    
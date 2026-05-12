"""
Unit tests for annotate_scan.py functions.

Run with: pytest test_annotate_scan.py -v
"""

import pytest
import numpy as np
import pandas as pd
from annotate_scan import (
    parse_formula, calculate_mass, add_formula_dicts, subtract_formula_dicts,
    formula_dict_to_string, get_isotope_distribution, cut_mmw_list,
    calcul_imp_mSigma, calcul_imp_formula
)


class TestFormulaParser:
    """Test formula parsing functions."""
    
    def test_parse_simple_formula(self):
        """Test parsing simple formulas."""
        result = parse_formula("C12H22O11")
        assert result == {'C': 12, 'H': 22, 'O': 11}
    
    def test_parse_formula_single_atoms(self):
        """Test parsing formulas with single atoms."""
        result = parse_formula("CHO")
        assert result == {'C': 1, 'H': 1, 'O': 1}
    
    def test_parse_formula_complex(self):
        """Test parsing complex oligonucleotide formula."""
        result = parse_formula("C189H238O119N66P18S4F8")
        assert result['C'] == 189
        assert result['H'] == 238
        assert result['O'] == 119
        assert result['N'] == 66
        assert result['P'] == 18
        assert result['S'] == 4
        assert result['F'] == 8
    
    def test_parse_empty_formula(self):
        """Test parsing empty/invalid formulas."""
        assert parse_formula("") == {}
        assert parse_formula("N/A") == {}
    
    def test_formula_dict_to_string(self):
        """Test converting formula dict back to string."""
        formula_dict = {'C': 12, 'H': 22, 'O': 11}
        result = formula_dict_to_string(formula_dict)
        assert result == "C12H22O11"
    
    def test_formula_dict_to_string_order(self):
        """Test that formula string has correct element order (C, H, rest)."""
        formula_dict = {'O': 5, 'H': 10, 'C': 6, 'N': 1}
        result = formula_dict_to_string(formula_dict)
        assert result.startswith("C6H10")


class TestMassCalculation:
    """Test mass calculation functions."""
    
    def test_calculate_monoisotopic_mass(self):
        """Test monoisotopic mass calculation."""
        formula = parse_formula("H2O")
        mass = calculate_mass(formula, 'mono')
        # H2O mono: 2*1.007825 + 15.994915 = 18.010565
        assert abs(mass - 18.010565) < 0.001
    
    def test_calculate_average_mass(self):
        """Test average mass calculation."""
        formula = parse_formula("H2O")
        mass = calculate_mass(formula, 'avg')
        # H2O avg: 2*1.00794 + 15.9994 = 18.01528
        assert abs(mass - 18.01528) < 0.001
    
    def test_empty_formula_mass(self):
        """Test mass of empty formula."""
        mass = calculate_mass({}, 'mono')
        assert mass == 0.0
    
    def test_add_formula_dicts(self):
        """Test formula addition."""
        f1 = {'C': 10, 'H': 20, 'O': 5}
        f2 = {'C': 2, 'H': 4, 'O': 1}
        result = add_formula_dicts(f1, f2)
        assert result == {'C': 12, 'H': 24, 'O': 6}
    
    def test_subtract_formula_dicts(self):
        """Test formula subtraction."""
        f1 = {'C': 10, 'H': 20, 'O': 5}
        f2 = {'C': 2, 'H': 4, 'O': 1}
        result = subtract_formula_dicts(f1, f2)
        assert result == {'C': 8, 'H': 16, 'O': 4}


class TestIsotopeDistribution:
    """Test isotope distribution functions."""
    
    def test_get_isotope_distribution_water(self):
        """Test isotope distribution for water."""
        formula = parse_formula("H2O")
        masses, intensities = get_isotope_distribution(formula, ntheo=3)
        
        # Should have isotope peaks
        assert len(masses) > 0
        assert len(intensities) > 0
        assert len(masses) == len(intensities)
        
        # Intensities should sum to approximately 1
        assert abs(np.sum(intensities) - 1.0) < 0.1
        
        # First mass should be close to monoisotopic mass
        expected_mono = calculate_mass(formula, 'mono')
        assert abs(masses[0] - expected_mono) < 0.1
    
    def test_get_isotope_distribution_empty(self):
        """Test isotope distribution for empty formula."""
        masses, intensities = get_isotope_distribution({}, ntheo=5)
        assert len(masses) == 0
        assert len(intensities) == 0


class TestClusteringFunctions:
    """Test molecular weight clustering."""
    
    def test_cut_mmw_list_single_cluster(self):
        """Test clustering with single cluster."""
        mw = np.array([1000.0, 1000.5, 1001.0])
        intensity = np.array([100, 150, 200])
        result = cut_mmw_list(mw, intensity, mw_window=10)
        
        assert np.all(result['id'] == 1)  # All in same cluster
        assert len(result['mw']) == 3
    
    def test_cut_mmw_list_multiple_clusters(self):
        """Test clustering with multiple clusters."""
        mw = np.array([1000.0, 1000.5, 1050.0, 1050.5])
        intensity = np.array([100, 150, 200, 250])
        result = cut_mmw_list(mw, intensity, mw_window=10)
        
        assert result['id'][0] == 1
        assert result['id'][1] == 1
        assert result['id'][2] == 2
        assert result['id'][3] == 2
    
    def test_cut_mmw_list_empty(self):
        """Test clustering with empty input."""
        result = cut_mmw_list(np.array([]), np.array([]), mw_window=10)
        assert len(result['id']) == 0
        assert len(result['mw']) == 0


class TestIsotopeMatching:
    """Test isotope pattern matching functions."""
    
    def test_calcul_imp_mSigma_valid_match(self):
        """Test isotope score calculation with valid data."""
        # Create synthetic experimental spectrum
        sp_deconvoluted = pd.DataFrame({
            'MW': [18.010, 18.012, 18.020],
            'I': [1.0, 0.3, 0.05]
        })
        
        # Create synthetic theoretical isotope pattern
        theo_masses = np.array([18.010, 18.012, 18.020])
        theo_intensities = np.array([1.0, 0.3, 0.05])
        
        result = calcul_imp_mSigma(sp_deconvoluted, theo_masses, theo_intensities,
                                   ntheo=3, max_mmw_ppm=10, baseline=100)
        
        assert result is not None
        assert 'score' in result
        assert 'oc_score' in result
        assert 'mono_mass' in result
        assert result['score'] >= 0
        assert 0 <= result['oc_score'] <= 1
    
    def test_calcul_imp_mSigma_empty_spectrum(self):
        """Test isotope score with empty spectrum."""
        sp_deconvoluted = pd.DataFrame({'MW': [], 'I': []})
        theo_masses = np.array([100.0])
        theo_intensities = np.array([1.0])
        
        result = calcul_imp_mSigma(sp_deconvoluted, theo_masses, theo_intensities,
                                   ntheo=1, max_mmw_ppm=10, baseline=100)
        
        assert result is None


class TestTransformationList:
    """Test transformation list expansion functions."""
    
    def test_calcul_imp_formula_single_transformation(self):
        """Test formula calculation for transformations."""
        trans_list = pd.DataFrame({
            'Plus_Formula': ['O'],
            'Minus_Formula': ['H2'],
            'CPD': ['Oxidation']
        })
        
        ifl_list = calcul_imp_formula("CH4", trans_list)
        
        assert len(ifl_list) == 1
        # CH4 + O - H2 = CH2O2 (but let's check counts)
        expected = {'C': 1, 'H': 2, 'O': 1}
        assert ifl_list[0] == expected
    
    def test_calcul_imp_formula_with_N_A(self):
        """Test formula calculation with N/A values."""
        trans_list = pd.DataFrame({
            'Plus_Formula': ['N/A'],
            'Minus_Formula': ['N/A'],
            'CPD': ['FLP']
        })
        
        ifl_list = calcul_imp_formula("C6H12O6", trans_list)
        
        assert len(ifl_list) == 1
        assert ifl_list[0] == {'C': 6, 'H': 12, 'O': 6}


class TestIntegration:
    """Integration tests combining multiple functions."""
    
    def test_full_pipeline_glucose(self):
        """Test full pipeline for glucose."""
        # Parse formula
        formula = parse_formula("C6H12O6")
        
        # Calculate masses
        mono_mass = calculate_mass(formula, 'mono')
        avg_mass = calculate_mass(formula, 'avg')
        
        assert 180.0 < mono_mass < 181.0
        assert 180.0 < avg_mass < 181.0
        
        # Get isotope distribution
        masses, intensities = get_isotope_distribution(formula, ntheo=5)
        
        assert len(masses) > 0
        assert abs(masses[0] - mono_mass) < 0.01
    
    def test_transformation_application(self):
        """Test applying transformations."""
        # Start with base formula
        base = parse_formula("C12H24O12")
        
        # Apply transformations
        trans1 = parse_formula("O")
        trans2 = parse_formula("H2")
        
        result = add_formula_dicts(base, trans1)
        result = subtract_formula_dicts(result, trans2)
        
        assert result['C'] == 12
        assert result['H'] == 22
        assert result['O'] == 13


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])

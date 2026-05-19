"""
Translation of R annotate_scan_targeted functions to Python.
Provides functions to annotate deconvoluted oligonucleotide spectra
based on transformation lists using isotope pattern matching.

Usage:
    from annotate_scan import annotate_scan_targeted
    result = annotate_scan_targeted(scan_processed_aggregated, 
                                    formula_flp='C192H239O117N73P18S4F8',
                                    cpd_flp='Demo A',
                                    transformation_list=trans_df,
                                    ntheo=12)

Dependencies: pandas, numpy, pyteomics
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple, List, Any
import re
import pyteomics.mass as mass

try:
    import IsoSpecPy
    HAS_ISOSPECPY = True
except ImportError:
    IsoSpecPy = None
    HAS_ISOSPECPY = False

def _isotopic_variants(formula_dict, npeaks=10):
    """
    Calculate isotopic variants for a given molecular formula.
    
    Args:
        formula_dict (dict): Dictionary with element symbols as keys and counts as values
        npeaks (int): Number of isotopic peaks to return
        
    Returns:
        list: List of (mass, abundance) tuples for the top npeaks isotopic variants
    """
    formula = formula_dict_to_string(formula_dict)
    if not formula:
        return []

    if HAS_ISOSPECPY:
        try:
            iso = IsoSpecPy.IsoThreshold(1e-12, formula=formula, absolute=False, get_confs=False)
            masses = np.array([iso.masses[i] for i in range(len(iso.masses))], dtype=float)
            intensities = np.array([iso.probs[i] for i in range(len(iso.probs))], dtype=float)

            if masses.size == 0:
                return []

            # Select top intensity peaks and sort by mass
            order = np.argsort(intensities)[::-1][:npeaks]
            masses = masses[order]
            intensities = intensities[order]
            sort_by_mass = np.argsort(masses)
            masses = masses[sort_by_mass]
            intensities = intensities[sort_by_mass]

            total_intensity = intensities.sum()
            if total_intensity > 0:
                intensities = intensities / total_intensity

            return list(zip(masses.tolist(), intensities.tolist()))
        except Exception:
            pass

    # Fallback to pyteomics isotopologues calculation
    formula = ''.join(f'{elem}{count}' for elem, count in formula_dict.items())
    iso_list = list(mass.isotopologues(formula, report_abundance=True, overall_threshold=1e-6))
    iso_list.sort(key=lambda x: x[1], reverse=True)
    top_iso = iso_list[:npeaks]
    result = []
    for comp, abund in top_iso:
        m = mass.calculate_mass(comp)
        result.append((m, abund))

    return result

# Atomic masses (for formula parsing and mass calculation)
ELEMENT_MASSES_MONO = {
    'H': 1.007825, 'C': 12.0, 'N': 14.003074, 'O': 15.994915,
    'P': 30.973762, 'S': 31.972072, 'F': 18.998403, 'Cl': 34.968853,
    'Br': 78.918336, 'I': 126.904477, 'Si': 27.976927, 'Sn': 119.902197,
    'B': 11.009305, 'Na': 22.989220, 'K': 38.963158, 'Fe': 55.934939
}

ELEMENT_MASSES_AVG = {
    'H': 1.00794, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
    'P': 30.973761, 'S': 32.065, 'F': 18.9984, 'Cl': 35.453,
    'Br': 79.904, 'I': 126.904, 'Si': 28.0855, 'Sn': 118.710,
    'B': 10.811, 'Na': 22.98977, 'K': 39.0983, 'Fe': 55.845
}


def parse_formula(formula_str: str) -> Dict[str, int]:
    """
    Parse elemental formula string (e.g., 'C12H22O11') into element counts.
    
    Parameters
    ----------
    formula_str : str
        Elemental formula string
        
    Returns
    -------
    dict
        Dictionary with element symbols as keys and counts as values
    """
    if pd.isna(formula_str):
        return {}

    if isinstance(formula_str, bytes):
        formula_str = formula_str.decode('utf-8', errors='ignore')

    if not isinstance(formula_str, str):
        formula_str = str(formula_str)

    formula_str = formula_str.strip()
    if not formula_str or formula_str == 'N/A':
        return {}
    
    element_counts = {}
    # Find all element-number pairs
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula_str)
    
    for element, count in matches:
        if element:  # Skip empty matches
            count = int(count) if count else 1
            element_counts[element] = element_counts.get(element, 0) + count
    
    return element_counts


def calculate_mass(formula_dict: Dict[str, int], mass_type: str = 'mono') -> float:
    """
    Calculate monoisotopic or average mass of a formula.
    
    Parameters
    ----------
    formula_dict : dict
        Dictionary with element symbols as keys and counts as values
    mass_type : str
        'mono' for monoisotopic mass, 'avg' for average mass
        
    Returns
    -------
    float
        Calculated mass
    """
    mass_table = ELEMENT_MASSES_MONO if mass_type == 'mono' else ELEMENT_MASSES_AVG
    total_mass = 0.0
    
    for element, count in formula_dict.items():
        if element in mass_table:
            total_mass += mass_table[element] * count
    
    return total_mass


def add_formula_dicts(formula1: Dict[str, int], formula2: Dict[str, int]) -> Dict[str, int]:
    """Add two formula dictionaries together."""
    result = formula1.copy()
    for element, count in formula2.items():
        result[element] = result.get(element, 0) + count
    return result


def subtract_formula_dicts(formula1: Dict[str, int], formula2: Dict[str, int]) -> Dict[str, int]:
    """Subtract formula2 from formula1."""
    result = formula1.copy()
    for element, count in formula2.items():
        result[element] = result.get(element, 0) - count
    return result


def formula_dict_to_string(formula_dict: Dict[str, int]) -> str:
    """Convert formula dictionary back to string format."""
    if not formula_dict:
        return ""
    
    # Sort elements: C, H, then rest alphabetically
    sorted_elements = []
    if 'C' in formula_dict:
        sorted_elements.append('C')
    if 'H' in formula_dict:
        sorted_elements.append('H')
    sorted_elements.extend(sorted([e for e in formula_dict.keys() if e not in ['C', 'H']]))
    
    result = ""
    for element in sorted_elements:
        count = formula_dict[element]
        if count > 0:
            if count == 1:
                result += element
            else:
                result += f"{element}{count}"
    
    return result


def get_isotope_distribution(formula_dict: Dict[str, int], ntheo: int = 12) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get theoretical isotope distribution using IsoSpecPy if available.
    
    Parameters
    ----------
    formula_dict : dict
        Dictionary with element symbols as keys and counts as values
    ntheo : int
        Number of theoretical isotope peaks
    
    Returns
    -------
    tuple
        (masses array, intensities array)
    """
    if not formula_dict:
        return np.array([]), np.array([])

    try:
        iso_variants = _isotopic_variants(formula_dict, npeaks=ntheo)
        if not iso_variants:
            return np.array([]), np.array([])

        masses = np.array([variant[0] for variant in iso_variants], dtype=float)
        intensities = np.array([variant[1] for variant in iso_variants], dtype=float)
        return masses, intensities
    except Exception as e:
        print(f"Error calculating isotope distribution: {e}")
        return np.array([]), np.array([])


def calcul_imp_mSigma(sp_deconvoluted: pd.DataFrame, theo_isotope_masses: np.ndarray,
                      theo_isotope_intensities: np.ndarray, ntheo: int, 
                      max_mmw_ppm: float, baseline: float) -> Optional[Dict]:
    """
    Calculate isotope pattern match score for a single compound.
    
    Parameters
    ----------
    sp_deconvoluted : pd.DataFrame
        Experimental spectrum with columns 'MW' and 'I'
    theo_isotope_masses : np.ndarray
        Theoretical isotope masses
    theo_isotope_intensities : np.ndarray
        Theoretical isotope intensities
    ntheo : int
        Number of theoretical peaks
    max_mmw_ppm : float
        Maximum allowed ppm error
    baseline : float
        Noise baseline level
        
    Returns
    -------
    dict or None
        Dictionary with matching scores and details, or None if no match
    """
    if len(theo_isotope_masses) == 0 or len(sp_deconvoluted) == 0:
        return None
    
    # Normalize theoretical distribution
    theo_deconvoluted = pd.DataFrame({
        'MW': theo_isotope_masses,
        'I': theo_isotope_intensities
    })
    theo_deconvoluted['I'] = theo_deconvoluted['I'] / theo_deconvoluted['I'].sum()
    
    tmw = theo_deconvoluted['MW'].iloc[0]  # Theoretical monoisotopic mass
    taw = (theo_deconvoluted['MW'] * theo_deconvoluted['I']).sum()  # Theoretical average mass
    
    sp_deconvoluted = sp_deconvoluted.copy()
    sp_deconvoluted.columns = ['MW', 'I']
    NT = len(theo_deconvoluted)
    
    # Initialize arrays for matching
    em_list = np.ones(NT)  # Error list
    res_list = np.full(NT, baseline / 2)  # Response list
    exp_list = np.zeros(NT)  # Experimental mass list
    
    abs_dev = 0.1  # Dalton tolerance
    
    # Match experimental to theoretical peaks
    for j in range(NT):
        errors = np.abs(sp_deconvoluted['MW'].values - theo_deconvoluted['MW'].iloc[j])
        valid = np.where(errors < abs_dev)[0]
        
        if len(valid) > 1:
            valid = np.array([valid[np.argmin(errors[valid])]])
        
        if len(valid) == 1:
            idx = valid[0]
            em_list[j] = errors[idx]
            res_list[j] = sp_deconvoluted['I'].iloc[idx]
            exp_list[j] = sp_deconvoluted['MW'].iloc[idx]
    
    # Calculate scores
    tm_list = theo_deconvoluted['MW'].values
    theo_list = theo_deconvoluted['I'].values
    idx = np.where(res_list > baseline)[0]
    
    if len(idx) == 0:
        return None
    
    # Overlap coefficient
    oc_score = len(idx) / NT
    
    res_list_norm = res_list / res_list.sum()
    theo_list_norm = theo_list / theo_list.sum()
    
    # Chi-square score
    chi_squa_score = np.sum(np.abs(res_list_norm[idx] - theo_list_norm[idx]) / theo_list_norm[idx]) / len(idx)
    
    # Calculate monoisotopic mass
    mono_mass = 0
    mono_ppm_dev = -1
    if len(idx) > 0 and idx[0] == 0:  # Monoisotopic peak detected
        mono_mass = exp_list[0]
        mono_ppm_dev = np.abs(em_list[0]) / tmw * 1000000
    
    # Calculate average mass
    avg_mass = np.sum(exp_list[idx] * res_list_norm[idx] / res_list_norm[idx].sum())
    taw_bis = np.sum(tm_list * theo_list_norm / theo_list_norm.sum())
    avg_mass_dev = np.abs(avg_mass - taw_bis)
    
    if avg_mass > 0:
        return {
            'score': chi_squa_score,
            'oc_score': oc_score,
            'mono_mass_ref': tmw,
            'mono_mass': mono_mass,
            'avg_mass_ref': taw_bis,
            'avg_mass': avg_mass,
            'mono_ppm_dev': mono_ppm_dev,
            'avg_mass_dev': avg_mass_dev,
            'exp_sp': np.column_stack([exp_list, res_list])
        }
    
    return None


def calcul_mix_mSigma(sp_deconvoluted: pd.DataFrame, theo_isotope1_masses: np.ndarray,
                      theo_isotope1_intensities: np.ndarray, theo_isotope2_masses: np.ndarray,
                      theo_isotope2_intensities: np.ndarray, coef1: float, coef2: float,
                      ntheo: int, max_mmw_ppm: float, baseline: float) -> Optional[Dict]:
    """
    Calculate isotope pattern match score for a mixture of two compounds.
    """
    if len(theo_isotope1_masses) == 0 or len(theo_isotope2_masses) == 0:
        return None
    
    # Normalize theoretical distributions
    theo_deconvoluted1 = pd.DataFrame({
        'MW': theo_isotope1_masses,
        'I': theo_isotope1_intensities
    })
    theo_deconvoluted1['I'] = theo_deconvoluted1['I'] / theo_deconvoluted1['I'].sum()
    
    theo_deconvoluted2 = pd.DataFrame({
        'MW': theo_isotope2_masses,
        'I': theo_isotope2_intensities
    })
    theo_deconvoluted2['I'] = theo_deconvoluted2['I'] / theo_deconvoluted2['I'].sum()
    
    tmw1 = theo_deconvoluted1['MW'].iloc[0]
    taw1 = (theo_deconvoluted1['MW'] * theo_deconvoluted1['I']).sum()
    
    tmw2 = theo_deconvoluted2['MW'].iloc[0]
    taw2 = (theo_deconvoluted2['MW'] * theo_deconvoluted2['I']).sum()
    
    # Mix the two distributions
    theo_deconvoluted1['I'] = theo_deconvoluted1['I'] * coef1
    theo_deconvoluted2['I'] = theo_deconvoluted2['I'] * coef2
    
    theo_mixed = pd.concat([theo_deconvoluted1, theo_deconvoluted2], ignore_index=True)
    theo_mixed = theo_mixed.sort_values('MW').reset_index(drop=True)
    
    sp_deconvoluted = sp_deconvoluted.copy()
    sp_deconvoluted.columns = ['MW', 'I']
    
    # Matching routine
    NT1 = len(theo_deconvoluted1)
    NT2 = len(theo_deconvoluted2)
    NT = len(theo_mixed)
    
    em_list1 = np.ones(NT1)
    em_list2 = np.ones(NT2)
    em_list = np.ones(NT)
    
    res_list1 = np.full(NT1, baseline / 2)
    res_list2 = np.full(NT2, baseline / 2)
    res_list = np.full(NT, baseline / 2)
    
    exp_list1 = np.zeros(NT1)
    exp_list2 = np.zeros(NT2)
    exp_list = np.zeros(NT)
    
    abs_dev = 0.2  # Higher tolerance for mixtures
    
    # Match experimental spectra to theoretical
    for j in range(NT1):
        errors = np.abs(sp_deconvoluted['MW'].values - theo_deconvoluted1['MW'].iloc[j])
        valid = np.where(errors < abs_dev)[0]
        if len(valid) > 1:
            valid = np.array([valid[np.argmin(errors[valid])]])
        if len(valid) == 1:
            idx = valid[0]
            em_list1[j] = errors[idx]
            res_list1[j] = sp_deconvoluted['I'].iloc[idx]
            exp_list1[j] = sp_deconvoluted['MW'].iloc[idx]
    
    for j in range(NT2):
        errors = np.abs(sp_deconvoluted['MW'].values - theo_deconvoluted2['MW'].iloc[j])
        valid = np.where(errors < abs_dev)[0]
        if len(valid) > 1:
            valid = np.array([valid[np.argmin(errors[valid])]])
        if len(valid) == 1:
            idx = valid[0]
            em_list2[j] = errors[idx]
            res_list2[j] = sp_deconvoluted['I'].iloc[idx]
            exp_list2[j] = sp_deconvoluted['MW'].iloc[idx]
    
    for j in range(NT):
        errors = np.abs(sp_deconvoluted['MW'].values - theo_mixed['MW'].iloc[j])
        valid = np.where(errors < abs_dev)[0]
        if len(valid) > 1:
            valid = np.array([valid[np.argmin(errors[valid])]])
        if len(valid) == 1:
            idx = valid[0]
            em_list[j] = errors[idx]
            res_list[j] = sp_deconvoluted['I'].iloc[idx]
            exp_list[j] = sp_deconvoluted['MW'].iloc[idx]
    
    # Calculate scores
    tm_list1 = theo_deconvoluted1['MW'].values
    theo_list1 = theo_deconvoluted1['I'].values
    tm_list2 = theo_deconvoluted2['MW'].values
    theo_list2 = theo_deconvoluted2['I'].values
    tm_list = theo_mixed['MW'].values
    theo_list = theo_mixed['I'].values
    
    idx1 = np.where(res_list1 > baseline)[0]
    idx2 = np.where(res_list2 > baseline)[0]
    idx = np.where(res_list > baseline)[0]
    
    oc_score1 = len(idx1) / NT1 if NT1 > 0 else 0
    oc_score2 = len(idx2) / NT2 if NT2 > 0 else 0
    oc_score = len(idx) / NT if NT > 0 else 0
    
    res_list1_norm = res_list1 / res_list1.sum()
    res_list2_norm = res_list2 / res_list2.sum()
    res_list_norm = res_list / res_list.sum()
    
    theo_list1_norm = theo_list1 / theo_list1.sum()
    theo_list2_norm = theo_list2 / theo_list2.sum()
    theo_list_norm = theo_list / theo_list.sum()
    
    chi_squa_score1 = np.sum(np.abs(res_list1_norm[idx1] - theo_list1_norm[idx1]) / theo_list1_norm[idx1]) / len(idx1) if len(idx1) > 0 else 0
    chi_squa_score2 = np.sum(np.abs(res_list2_norm[idx2] - theo_list2_norm[idx2]) / theo_list2_norm[idx2]) / len(idx2) if len(idx2) > 0 else 0
    chi_squa_score = np.sum(np.abs(res_list_norm[idx] - theo_list_norm[idx]) / theo_list_norm[idx]) / len(idx) if len(idx) > 0 else 0
    
    # Calculate masses
    mono_mass1 = 0
    mono_ppm_dev1 = -1
    if len(idx1) > 0 and idx1[0] == 0:
        mono_mass1 = exp_list1[0]
        mono_ppm_dev1 = np.abs(em_list1[0]) / tmw1 * 1000000 if tmw1 > 0 else -1
    
    mono_mass2 = 0
    mono_ppm_dev2 = -1
    if len(idx2) > 0 and idx2[0] == 0:
        mono_mass2 = exp_list2[0]
        mono_ppm_dev2 = np.abs(em_list2[0]) / tmw2 * 1000000 if tmw2 > 0 else -1
    
    avg_mass1 = np.sum(exp_list1[idx1] * theo_list1_norm[idx1] / theo_list1_norm[idx1].sum()) if len(idx1) > 0 else 0
    avg_mass2 = np.sum(exp_list2[idx2] * theo_list2_norm[idx2] / theo_list2_norm[idx2].sum()) if len(idx2) > 0 else 0
    avg_mass = np.sum(exp_list[idx] * theo_list_norm[idx] / theo_list_norm[idx].sum()) if len(idx) > 0 else 0
    
    taw_bis1 = np.sum(tm_list1 * theo_list1_norm / theo_list1_norm.sum())
    taw_bis2 = np.sum(tm_list2 * theo_list2_norm / theo_list2_norm.sum())
    taw_bis = np.sum(tm_list * theo_list_norm / theo_list_norm.sum())
    
    avg_mass_dev1 = np.abs(avg_mass1 - taw_bis1)
    avg_mass_dev2 = np.abs(avg_mass2 - taw_bis2)
    avg_mass_dev = np.abs(avg_mass - taw_bis)
    
    return {
        'score1': chi_squa_score1, 'score2': chi_squa_score2, 'score': chi_squa_score,
        'oc_score1': oc_score1, 'oc_score2': oc_score2, 'oc_score': oc_score,
        'mono_mass_ref1': tmw1, 'avg_mass_ref1': taw_bis1,
        'mono_mass1': mono_mass1, 'avg_mass1': avg_mass1,
        'mono_mass_ref2': tmw2, 'avg_mass_ref2': taw_bis2,
        'mono_mass2': mono_mass2, 'avg_mass2': avg_mass2,
        'avg_mass_ref': taw_bis, 'avg_mass': avg_mass,
        'avg_mass_dev': avg_mass_dev,
        'mono_ppm_dev1': mono_ppm_dev1, 'avg_mass_dev1': avg_mass_dev1,
        'mono_ppm_dev2': mono_ppm_dev2, 'avg_mass_dev2': avg_mass_dev2,
        'exp_sp1': np.column_stack([exp_list1, res_list1]),
        'exp_sp2': np.column_stack([exp_list2, res_list2]),
        'exp_sp': np.column_stack([exp_list, res_list])
    }


def cut_mmw_list(mwlist: np.ndarray, intlist: np.ndarray, mw_window: float) -> Dict:
    """
    Cluster molecular weights into features based on mass window.
    
    Parameters
    ----------
    mwlist : np.ndarray
        Array of molecular weights
    intlist : np.ndarray
        Array of intensities
    mw_window : float
        Mass window for clustering
        
    Returns
    -------
    dict
        Dictionary with 'id' (feature IDs) and 'mw' (average mass per feature)
    """
    mwlist = np.asarray(mwlist)
    intlist = np.asarray(intlist)
    N = len(mwlist)
    
    if N == 0:
        return {'id': np.array([], dtype=int), 'mw': np.array([])}
    
    f = 1
    mw_feature = np.zeros(N, dtype=int)
    mw_avg = np.zeros(N, dtype=float)
    t0 = 0
    
    for k in range(1, N):
        ttt = np.arange(t0, k)
        min_mw = np.min(mwlist[ttt])
        idx_best = np.argmax(intlist[ttt])
        best_mw = mwlist[ttt[idx_best]]
        max_mw = np.max(mwlist[ttt])
        
        if (mwlist[k] - min_mw > mw_window) or (mwlist[k] - max_mw > 0.2):
            mw_feature[t0:k] = f
            mw_avg[t0:k] = np.round(best_mw, 4)
            f += 1
            t0 = k
    
    ttt = np.arange(t0, N)
    mw_feature[ttt] = f
    idx_best = ttt[np.argmax(intlist[ttt])] if len(ttt) > 0 else 0
    mw_avg[ttt] = np.round(mwlist[idx_best], 4) if len(ttt) > 0 else 0.0
    
    return {'id': mw_feature, 'mw': mw_avg}


def calcul_imp_formula(formula_flp: str, transformation_list: pd.DataFrame) -> List[Dict[str, int]]:
    """
    Calculate impurity formulas by applying transformation list.
    
    Parameters
    ----------
    formula_flp : str
        Formula of full length product
    transformation_list : pd.DataFrame
        DataFrame with Plus_Formula and Minus_Formula columns
        
    Returns
    -------
    list
        List of formula dictionaries for each transformation
    """
    ifl_list = []
    
    if formula_flp and formula_flp != 'N/A':
        rfl = parse_formula(formula_flp)
        
        for idx, row in transformation_list.iterrows():
            plus_formula = row.get('Plus_Formula', 'N/A')
            minus_formula = row.get('Minus_Formula', 'N/A')
            
            if pd.isna(plus_formula) or plus_formula == 'N/A':
                pfl = {}
            else:
                pfl = parse_formula(plus_formula)

            if pd.isna(minus_formula) or minus_formula == 'N/A':
                mfl = {}
            else:
                mfl = parse_formula(minus_formula)
            
            ifl = add_formula_dicts(rfl, pfl)
            ifl = subtract_formula_dicts(ifl, mfl)
            ifl_list.append(ifl)
    else:
        # Use formula from transformation list directly
        for idx, row in transformation_list.iterrows():
            formula = row.get('FORMULA', '')
            ifl = parse_formula(formula)
            ifl_list.append(ifl)
    
    return ifl_list


def expand_transformation_list(formula_flp: str, transformation_list: pd.DataFrame) -> Tuple[List, pd.DataFrame]:
    """
    Expand transformation list with formulas and masses.
    
    Parameters
    ----------
    formula_flp : str
        Formula of full length product
    transformation_list : pd.DataFrame
        Base transformation list
        
    Returns
    -------
    tuple
        (IFL list, expanded transformation_list)
    """
    # Calculate impurity formulas
    IFL = calcul_imp_formula(formula_flp, transformation_list)
    
    # Calculate masses
    amw_flp = calculate_mass(parse_formula(formula_flp), 'avg') if formula_flp and formula_flp != 'N/A' else 0
    mmw_flp = calculate_mass(parse_formula(formula_flp), 'mono') if formula_flp and formula_flp != 'N/A' else 0
    
    # Add formula and mass information to transformation list
    transformation_list = transformation_list.copy()
    transformation_list['FORMULA'] = [formula_dict_to_string(ifl) for ifl in IFL]
    transformation_list['AVG.MW'] = [calculate_mass(ifl, 'avg') for ifl in IFL]
    transformation_list['MONO.MW'] = [calculate_mass(ifl, 'mono') for ifl in IFL]
    transformation_list['Delta.AVG.MW'] = transformation_list['AVG.MW'] - amw_flp
    transformation_list['Delta.MONO.MW'] = transformation_list['MONO.MW'] - mmw_flp
    
    return IFL, transformation_list


def annotate_envelop(envelop: pd.DataFrame, ref_trans: pd.DataFrame, IFL: List[Dict[str, int]],
                     ntheo: int = 12, baseline: float = 1000, min_overlap: float = 0.6,
                     max_mmw_ppm: float = 10) -> Dict:
    """
    Annotate an isotope envelope by matching to reference transformations.
    
    Parameters
    ----------
    envelop : pd.DataFrame
        Envelope peaks with columns 'MW' and 'Response'
    ref_trans : pd.DataFrame
        Reference transformation list
    IFL : list
        List of impurity formulas
    ntheo : int
        Number of theoretical isotope peaks
    baseline : float
        Noise baseline level
    min_overlap : float
        Minimum overlap ratio for matching
    max_mmw_ppm : float
        Maximum ppm error allowed
        
    Returns
    -------
    dict
        Dictionary with 'envelop.annotated' and 'features.annotated'
    """
    envelop_annotated = envelop.copy()
    sd1 = envelop.copy()
    
    # Initialize annotation columns with object dtype so string labels can be assigned
    sd1['CPD'] = pd.Series([None] * len(sd1), index=sd1.index, dtype='object')
    sd1['FORMULA'] = pd.Series([None] * len(sd1), index=sd1.index, dtype='object')
    sd1['SCORE'] = pd.Series([None] * len(sd1), index=sd1.index, dtype='object')
    sd1['OC_SCORE'] = pd.Series([None] * len(sd1), index=sd1.index, dtype='object')
    sd1['CPD1'] = pd.Series([None] * len(sd1), index=sd1.index, dtype='object')
    
    # Calculate envelope center and range
    tmp_scan = sd1[['MW', 'Response']].copy()
    tmp_scan.columns = ['MW', 'I']
    
    total_intensity = tmp_scan['I'].sum()
    if total_intensity > 0:
        avg_envelop = (tmp_scan['MW'] * tmp_scan['I']).sum() / total_intensity
    else:
        avg_envelop = tmp_scan['MW'].mean()
    
    max_envelop = tmp_scan['MW'].max()
    min_envelop = tmp_scan['MW'].min()
    
    # Find candidate formulas within mass window
    valid = []
    avg_dev = []
    formula_valid = []
    kt = []
    
    for k in range(len(ref_trans)):
        dev_mass = abs(ref_trans['AVG.MW'].iloc[k] - avg_envelop)
        
        if dev_mass <= 10:
            # Get theoretical isotope pattern
            theo_masses, theo_intensities = get_isotope_distribution(IFL[k], ntheo)
            
            if len(theo_masses) == 0:
                continue
            
            # Check nominal mass overlap
            matched_nm = np.intersect1d(np.round(theo_masses), np.round(tmp_scan['MW'].values))
            
            if (ntheo > 6 and len(matched_nm) >= ntheo * min_overlap) or \
               (ntheo <= 6 and len(matched_nm) >= 3):
                kt.append(len(matched_nm))
                valid.append(k)
                avg_dev.append(dev_mass)
                formula_valid.append(ref_trans['FORMULA'].iloc[k])
    
    # Remove duplicate formulas
    if len(valid) > 1:
        unique_formulas, indices = np.unique(formula_valid, return_index=True)
        if len(indices) < len(valid):
            valid = [valid[i] for i in indices]
            avg_dev = [avg_dev[i] for i in indices]
            kt = [kt[i] for i in indices]
    
    coef1 = 1
    coef2 = 0
    mSigma = None
    
    # Single possibility
    if len(valid) == 1:
        ifl = IFL[valid[0]]
        theo_masses, theo_intensities = get_isotope_distribution(ifl, ntheo)
        mSigma = calcul_imp_mSigma(tmp_scan, theo_masses, theo_intensities, 
                                   ntheo, max_mmw_ppm, baseline)
        coef1 = 1
        coef2 = 0
    
    # Two or more possibilities
    elif len(valid) >= 2:
        # Rank by mass deviation and number of matched peaks
        tmp_r = np.argsort(avg_dev) + np.argsort([-k for k in kt])
        best_indices = np.argsort(tmp_r)[:2]
        valid = [valid[i] for i in best_indices]
        
        # Try to distinguish between the two
        ifl1 = IFL[valid[0]]
        ifl2 = IFL[valid[1]]
        
        theo_masses1, theo_intensities1 = get_isotope_distribution(ifl1, ntheo)
        theo_masses2, theo_intensities2 = get_isotope_distribution(ifl2, ntheo)
        
        # Try mixture deconvolution
        # Simple approach: assume equal coefficients first
        mSigma_single1 = calcul_imp_mSigma(tmp_scan, theo_masses1, theo_intensities1,
                                           ntheo, max_mmw_ppm, baseline)
        mSigma_single2 = calcul_imp_mSigma(tmp_scan, theo_masses2, theo_intensities2,
                                           ntheo, max_mmw_ppm, baseline)
        
        # Choose the better single match
        if mSigma_single1 and mSigma_single2:
            if mSigma_single1['score'] <= mSigma_single2['score']:
                mSigma = mSigma_single1
                coef1 = 1
                coef2 = 0
            else:
                mSigma = mSigma_single2
                coef1 = 0
                coef2 = 1
                valid = valid[1:2]
        elif mSigma_single1:
            mSigma = mSigma_single1
            coef1 = 1
            coef2 = 0
        elif mSigma_single2:
            mSigma = mSigma_single2
            coef1 = 0
            coef2 = 1
            valid = valid[1:2]
    
    # Add annotations to envelop
    if mSigma is not None and not np.isnan(mSigma.get('score', np.nan)):
        exp_sp = mSigma.get('exp_sp', np.array([]))
        if len(exp_sp) > 0:
            for pos in range(len(sd1)):
                row_index = sd1.index[pos]
                mw = sd1['MW'].iloc[pos]
                # Find matching theoretical mass
                matches = np.where(np.isclose(exp_sp[:, 0], mw, atol=0.1))[0]
                if len(matches) > 0:
                    sd1.at[row_index, 'CPD'] = ref_trans['CPD'].iloc[valid[0]]
                    sd1.at[row_index, 'FORMULA'] = ref_trans['FORMULA'].iloc[valid[0]]
                    sd1.at[row_index, 'SCORE'] = round(mSigma['score'], 2)
                    sd1.at[row_index, 'OC_SCORE'] = round(mSigma['oc_score'], 2)
                    sd1.at[row_index, 'CPD1'] = ref_trans['CPD'].iloc[valid[0]]
    
    # Clean up NaN values
    for col in ['CPD', 'FORMULA', 'SCORE', 'OC_SCORE']:
        sd1[col] = sd1[col].astype(str).str.replace('nan', '', regex=False)
        sd1[col] = sd1[col].astype(str).str.replace('N/A:', '', regex=False)
    
    sd1['COEF'] = f"{coef1:.2f}:{coef2:.2f}"
    
    # Create features summary
    features = []
    if mSigma is not None and not np.isnan(mSigma.get('score', np.nan)):
        feature_row = {
            'FEATURE': ref_trans['CPD'].iloc[valid[0]],
            'FORMULA': ref_trans['FORMULA'].iloc[valid[0]],
            'THEO_MMW': round(mSigma.get('mono_mass_ref', 0), 4),
            'THEO_AMW': round(mSigma.get('avg_mass_ref', 0), 4),
            'EXP_MMW': round(mSigma.get('mono_mass', 0), 4),
            'EXP_AMW': round(mSigma.get('avg_mass', 0), 4),
            'EXP_MMW_PPM': round(mSigma.get('mono_ppm_dev', -1), 2),
            'EXP_AMW_DEV': round(mSigma.get('avg_mass_dev', 0), 2),
            'SCORE': round(mSigma['score'], 2),
            'OC': round(mSigma['oc_score'], 2),
            'RESPONSE': round(sd1['Response'].sum(), 0),
            'Envelop': sd1['Envelop'].iloc[0] if 'Envelop' in sd1.columns else 1
        }
        features.append(feature_row)
    
    features_df = pd.DataFrame(features) if features else pd.DataFrame()
    
    return {
        'envelop.annotated': sd1,
        'features.annotated': features_df
    }


def annotate_scan_targeted(scan_processed_aggregated: Optional[pd.DataFrame] = None,
                          formula_flp: str = "C192H239O117N73P18S4F8",
                          cpd_flp: str = "Demo A",
                          transformation_list: Optional[pd.DataFrame] = None,
                          mdb: Optional[pd.DataFrame] = None,
                          ntheo: int = 12,
                          min_overlap: float = 0.6,
                          max_msigma: float = 3,
                          max_mmw_ppm: float = 10,
                          baseline: float = 1000) -> Dict[str, Optional[pd.DataFrame]]:
    """
    Annotate deconvoluted oligonucleotide spectra based on transformation list.
    
    Parameters
    ----------
    scan_processed_aggregated : pd.DataFrame
        Deconvoluted spectrum from process_scan output with columns
        'MW', 'Response', 'Envelop', etc.
    formula_flp : str
        Neutral elemental formula of main compound
    cpd_flp : str
        Name of main compound
    transformation_list : pd.DataFrame
        DataFrame with columns: ID, CPD, Plus_Formula, Minus_Formula,
        Delta.AVG.MW, Delta.MONO.MW
    mdb : pd.DataFrame
        Alternative: direct list of compounds to search without transformations
    ntheo : int
        Number of theoretical isotope peaks
    min_overlap : float
        Minimum overlap ratio (0-1) for isotope pattern matching
    max_msigma : float
        Maximum chi-square score allowed
    max_mmw_ppm : float
        Maximum ppm error allowed for mass matching
    baseline : float
        Noise baseline level
        
    Returns
    -------
    dict
        Dictionary with 'scan' (annotated spectrum) and 'feature' (summary table)
    """
    if scan_processed_aggregated is None or len(scan_processed_aggregated) == 0:
        return {'scan': None, 'feature': None}
    
    scan_annotated = None
    features_annotated = None
    
    # Prepare transformation list
    if transformation_list is None and mdb is None:
        return {'scan': None, 'feature': None}
    
    if transformation_list is not None and len(transformation_list) > 0:
        # Expand transformation list
        IFL, transformation_list = expand_transformation_list(formula_flp, transformation_list)
    elif mdb is not None:
        # Use molecular database directly
        transformation_list = mdb
        IFL = [parse_formula(row.get('FORMULA', '')) for _, row in mdb.iterrows()]
    else:
        return {'scan': None, 'feature': None}
    
    if len(transformation_list) == 0:
        return {'scan': None, 'feature': None}
    
    # Annotate each envelope
    if 'Envelop' in scan_processed_aggregated.columns:
        EEE = int(scan_processed_aggregated['Envelop'].max())
    else:
        EEE = 1
    
    results_list = []
    features_list = []
    
    for i in range(1, EEE + 1):
        inds = np.where(scan_processed_aggregated['Envelop'].values == i)[0]
        
        if len(inds) >= 2:  # At least 2 peaks in envelope
            envelop = scan_processed_aggregated.iloc[inds].copy()
            results = annotate_envelop(envelop, transformation_list, IFL,
                                       ntheo=ntheo, baseline=baseline,
                                       min_overlap=min_overlap, max_mmw_ppm=max_mmw_ppm)
            
            if results['envelop.annotated'] is not None:
                results_list.append(results['envelop.annotated'])
            if results['features.annotated'] is not None and len(results['features.annotated']) > 0:
                features_list.append(results['features.annotated'])
    
    if results_list:
        scan_annotated = pd.concat(results_list, ignore_index=True)
    
    if features_list:
        features_annotated = pd.concat(features_list, ignore_index=True)
        
        # Post-process: Keep best scoring features
        if len(features_annotated) > 0:
            features_annotated = features_annotated[features_annotated['RESPONSE'] > 0].copy()
            features_annotated = features_annotated.sort_values('SCORE').reset_index(drop=True)
            
            # Deduplicate by feature
            tmp_features = features_annotated['FEATURE'].unique()
            final_features = []
            
            for tf in tmp_features:
                valid = features_annotated[features_annotated['FEATURE'] == tf]
                if len(valid) > 0:
                    new_feature = valid.iloc[0].copy()
                    new_feature['RESPONSE'] = valid['RESPONSE'].sum()
                    final_features.append(new_feature)
            
            features_annotated = pd.DataFrame(final_features)
            
            # Filter by score and overlap
            features_annotated = features_annotated[
                (features_annotated['SCORE'] <= max_msigma) &
                (features_annotated['OC'] >= min_overlap)
            ].copy()
    
    return {
        'scan': scan_annotated,
        'feature': features_annotated
    }

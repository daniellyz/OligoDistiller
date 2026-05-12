# Python Translation of R `annotate_scan_targeted` - Documentation

## Overview

This is a complete Python translation of the R `annotate_scan_targeted` function from the OligoDistiller package. It annotates deconvoluted oligonucleotide mass spectra by matching experimental isotope patterns against theoretical patterns from a transformation list.

## Key Changes from R to Python

### 1. **BRAIN Package Replacement**
- **R**: Uses `BRAIN::useBRAIN()`, `calculateAverageMass()`, `calculateMonoisotopicMass()`
- **Python**: Uses `pyteomics` package via `mass.isotopologues()` function
  - Installation: `pip install brainpy`
  - Reference: https://github.com/mobiusklein/brainpy

### 2. **Formula Parsing**
- **R**: Uses regex-based `ListFormula1()` function
- **Python**: Implemented as `parse_formula()` function using standard regex

### 3. **Data Structures**
- **R**: Lists and data frames
- **Python**: Dictionaries and pandas DataFrames

### 4. **String Manipulation**
- **R**: `stringr` functions (str_trim, str_remove)
- **Python**: Built-in string methods

## Function Reference

### Main Function

#### `annotate_scan_targeted()`

Annotates a deconvoluted spectrum against a transformation list.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scan_processed_aggregated` | DataFrame | None | Output from `process_scan()` with aggregated peaks |
| `formula_flp` | str | "C192H239O117N73P18S4F8" | Elemental formula of full-length product |
| `cpd_flp` | str | "Demo A" | Name/ID of main compound |
| `transformation_list` | DataFrame | None | Transformation list with columns: ID, CPD, Plus_Formula, Minus_Formula, Delta.AVG.MW, Delta.MONO.MW |
| `mdb` | DataFrame | None | Molecular database (alternative to transformation_list) |
| `ntheo` | int | 12 | Number of theoretical isotope peaks |
| `min_overlap` | float | 0.6 | Minimum isotope pattern overlap (0-1) |
| `max_msigma` | float | 3 | Maximum chi-square score allowed |
| `max_mmw_ppm` | float | 10 | Maximum ppm error for mass matching |
| `baseline` | float | 1000 | Noise baseline level |

**Returns:**

Dictionary with keys:
- `'scan'`: Annotated spectrum DataFrame with columns: MW, Response, CPD, FORMULA, SCORE, OC_SCORE, COEF
- `'feature'`: Feature summary DataFrame with columns: FEATURE, FORMULA, THEO_MMW, THEO_AMW, EXP_MMW, EXP_AMW, EXP_MMW_PPM, EXP_AMW_DEV, SCORE, OC, RESPONSE, Envelop

**Example:**

```python
from annotate_scan import annotate_scan_targeted
import pandas as pd

# Load transformation list
trans_list = pd.read_csv('transformation_list.txt', sep='\t')

# Annotate spectrum
result = annotate_scan_targeted(
    scan_processed_aggregated=aggregated_spectrum,
    formula_flp="C189H238O119N66P18S4F8",
    cpd_flp="Demo A",
    transformation_list=trans_list,
    ntheo=10
)

print(result['scan'])   # Annotated spectrum
print(result['feature'])  # Summary features
```

### Helper Functions

#### `parse_formula(formula_str: str) -> Dict[str, int]`

Parses elemental formula string into element-count dictionary.

```python
parse_formula("C12H22O11")  # {'C': 12, 'H': 22, 'O': 11}
```

#### `calculate_mass(formula_dict: Dict[str, int], mass_type: str = 'mono') -> float`

Calculates monoisotopic or average mass of a formula.

```python
formula = parse_formula("C12H22O11")
mono_mass = calculate_mass(formula, 'mono')  # Monoisotopic mass
avg_mass = calculate_mass(formula, 'avg')    # Average mass
```

#### `get_isotope_distribution(formula_dict: Dict[str, int], ntheo: int = 12) -> Tuple[np.ndarray, np.ndarray]`

Gets theoretical isotope pattern using brainpy.

```python
masses, intensities = get_isotope_distribution(formula, ntheo=10)
```

#### `calcul_imp_mSigma(sp_deconvoluted, theo_isotope_masses, theo_isotope_intensities, ntheo, max_mmw_ppm, baseline)`

Calculates isotope pattern match score for a single compound.

**Returns**: Dictionary with keys:
- `score`: Chi-square score (lower is better)
- `oc_score`: Overlap coefficient
- `mono_mass_ref`: Theoretical monoisotopic mass
- `mono_mass`: Experimental monoisotopic mass
- `avg_mass_ref`: Theoretical average mass
- `avg_mass`: Experimental average mass
- `mono_ppm_dev`: PPM deviation for monoisotopic peak
- `avg_mass_dev`: Absolute deviation for average mass
- `exp_sp`: Matched experimental spectrum

#### `calcul_mix_mSigma(...)`

Calculates isotope pattern match for a mixture of two compounds.

**Returns**: Extended dictionary with separate scores and masses for each component (score1, score2, mono_mass1, mono_mass2, etc.)

#### `cut_mmw_list(mwlist, intlist, mw_window) -> Dict`

Clusters molecular weights into features based on mass window.

```python
clusters = cut_mmw_list(
    mwlist=np.array([1000.5, 1001.2, 1050.0]),
    intlist=np.array([100, 150, 200]),
    mw_window=10
)
# Returns: {'id': [1, 1, 2], 'mw': [1000.5, 1000.5, 1050.0]}
```

#### `calcul_imp_formula(formula_flp, transformation_list) -> List[Dict[str, int]]`

Calculates impurity formulas by applying transformations to full-length product.

```python
ifl_list = calcul_imp_formula("C189H238O119N66P18S4F8", transformation_list)
```

#### `expand_transformation_list(formula_flp, transformation_list) -> Tuple[List, DataFrame]`

Expands transformation list with calculated formulas and masses.

```python
ifl, expanded_trans = expand_transformation_list(formula_flp, trans_list)
```

#### `annotate_envelop(envelop, ref_trans, IFL, ntheo, baseline, min_overlap, max_mmw_ppm)`

Annotates a single isotope envelope.

## Transformation List Format

CSV/TSV file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| ID | int | Unique identifier |
| CPD | str | Compound name (e.g., "Phosphorothioate", "Deamination") |
| Plus_Formula | str | Formula to add (use "N/A" for none) |
| Minus_Formula | str | Formula to subtract (use "N/A" for none) |
| Delta.AVG.MW | float | Average mass difference |
| Delta.MONO.MW | float | Monoisotopic mass difference |

**Example:**

```
ID	CPD	Plus_Formula	Minus_Formula	Delta.AVG.MW	Delta.MONO.MW
1	FLP	N/A	N/A	0	0
4	Phosphorothioate	S	O	16.0667	15.9772
24	Dehydration	N/A	H2O	-18.0153	-18.0106
```

Download example: https://raw.githubusercontent.com/daniellyz/MESSAR/refs/heads/master/MESSAR_WEBSERVER_DEMO/Transformation_list_jennifer_shortened.txt

## Input/Output from process_scan

This function is designed to work with output from `process_scan()`:

### Input: `scan_processed_aggregated`

Required columns:
- **MW**: Molecular weight (float)
- **Response**: Intensity (float)
- **Envelop**: Envelope ID (int)

Optional columns:
- **Mass**, **z**, etc. (preserved in output)

### Output

#### `scan` (Annotated Spectrum)
Columns added/modified:
- **CPD**: Matched compound name
- **FORMULA**: Matched elemental formula
- **SCORE**: Chi-square match score
- **OC_SCORE**: Overlap coefficient
- **CPD1**: Cleaned compound name
- **COEF**: Mixture coefficients (coef1:coef2)

#### `feature` (Feature Summary)
- **FEATURE**: Feature name/ID
- **FORMULA**: Elemental formula
- **THEO_MMW**: Theoretical monoisotopic mass
- **THEO_AMW**: Theoretical average mass
- **EXP_MMW**: Experimental monoisotopic mass
- **EXP_AMW**: Experimental average mass
- **EXP_MMW_PPM**: PPM error for monoisotopic mass
- **EXP_AMW_DEV**: Absolute error for average mass
- **SCORE**: Chi-square score
- **OC**: Overlap coefficient
- **RESPONSE**: Total response intensity
- **Envelop**: Envelope ID

## Algorithm Overview

1. **Transformation Expansion**: Calculate elemental formulas for all transformations by:
   - Base formula (FLP) + Plus_Formula - Minus_Formula
   - Calculate monoisotopic and average masses

2. **Envelope Clustering**: Group peaks into envelopes if they have ≥2 peaks

3. **Formula Matching**: For each envelope, find candidate formulas within 10 Da of measured average mass

4. **Isotope Distribution**: Calculate theoretical isotope patterns for candidates using brainpy

5. **Pattern Matching**: 
   - For single candidate: Direct match
   - For multiple candidates: Evaluate both single and mixture scenarios

6. **Scoring**: Calculate chi-square score between experimental and theoretical patterns

7. **Filtering**: Keep only matches with score ≤ max_msigma and overlap ≥ min_overlap

## Important Parameters

### `ntheo` - Number of Isotope Peaks
- Higher values capture more isotope details but may be noisy
- Typical values: 6-12
- Adjust based on expected oligonucleotide mass and instrument resolution

### `min_overlap` - Isotope Pattern Overlap
- Fraction of theoretical peaks that must match
- Range: 0-1
- For tight matching: 0.6-0.8
- For loose matching: 0.3-0.5

### `max_msigma` - Chi-Square Score Threshold
- Lower = stricter matching (fewer false positives)
- Higher = more permissive (may catch noisy data)
- Typical: 3-5
- For noisy data: 10-20

### `max_mmw_ppm` - Mass Accuracy
- Maximum allowed ppm error
- Typically: 5-10 ppm for high-resolution MS

### `baseline` - Noise Level
- For MS1: Usually 1000
- For MS/MS: Usually 50-100
- Adjust to your instrument noise characteristics

## Performance Considerations

1. **Memory**: Linear with spectrum size
2. **Speed**: O(n × m) where n = spectrum peaks, m = transformations
3. **For large datasets**: Process in batches or use vectorized operations

## Error Handling

- Missing transformation list: Returns None for both outputs
- Empty spectrum: Returns None
- Invalid formulas: Skipped with warning
- Failed isotope calculation: Skipped with warning

## Dependencies

- **pandas**: Data manipulation
- **numpy**: Array operations
- **brainpy**: Isotope pattern calculation
  - Install: `pip install brainpy`

## References

- Original R code: https://github.com/daniellyz/OligoDistiller
- BRAIN algorithm paper: Valkenborg, D., et al. (2012)
- brainpy documentation: https://github.com/mobiusklein/brainpy

## Differences from R Implementation

1. **Automatic variable expansion**: Python version automatically expands transformation lists (no need for separate expand_transformation_list call)

2. **Simplified mixture detection**: Current version uses simpler logic for detecting 2-compound mixtures (compares single matches rather than full optimization)

3. **Data types**: Uses NumPy arrays and Pandas DataFrames instead of R vectors/lists

4. **String handling**: Cleaner regex-based formula parsing

5. **No error suppression**: Python version provides explicit error messages (R version had `options(error = expression(NULL))`)

## Future Improvements

1. Implement full mixture optimization using scipy.optimize
2. Add parallel processing for large transformation lists
3. Support multiple input formats (mzML, mzXML)
4. Add visualization functions
5. Implement caching for isotope calculations

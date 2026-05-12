# Quick Reference Guide for annotate_scan_targeted

## Installation

```bash
pip install -r requirements.txt
```

This installs:
- pandas, numpy for data handling
- brainpy for isotope pattern calculations
- pytest for testing

## 5-Minute Quick Start

### 1. Import and Prepare Data

```python
import pandas as pd
from process_scan import process_scan
from annotate_scan import annotate_scan_targeted

# Load transformation list
trans_list = pd.read_csv('transformation_list.txt', sep='\t')

# Or create manually:
trans_list = pd.DataFrame({
    'ID': [1, 2, 3],
    'CPD': ['FLP', 'Phosphorothioate', 'Dehydration'],
    'Plus_Formula': ['N/A', 'S', 'N/A'],
    'Minus_Formula': ['N/A', 'O', 'H2O'],
    'Delta.AVG.MW': [0, 16.0667, -18.0153],
    'Delta.MONO.MW': [0, 15.9772, -18.0106]
})
```

### 2. Process Raw Spectrum

```python
# Your raw MS data as (mass, intensity) pairs
raw_data = [(500.5, 1000), (501.2, 1500), ...]

result = process_scan(
    test_scan=raw_data,
    polarity='Negative',
    baseline=1000,
    min_mw=4000,
    max_mw=12000
)

aggregated = result['scan_processed_aggregated']
```

### 3. Annotate Spectrum

```python
annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    formula_flp="C189H238O119N66P18S4F8",
    cpd_flp="MyOligo",
    transformation_list=trans_list,
    ntheo=10,
    min_overlap=0.6,
    max_msigma=5,
    max_mmw_ppm=10,
    baseline=1000
)

# Extract results
spectrum_with_annotations = annotated['scan']
feature_summary = annotated['feature']

# Save results
spectrum_with_annotations.to_csv('annotated_spectrum.csv')
feature_summary.to_csv('annotated_features.csv')
```

### 4. View Results

```python
# Annotated peaks
print(spectrum_with_annotations[['MW', 'Response', 'CPD', 'FORMULA', 'SCORE']])

# Feature summary
print(feature_summary[['FEATURE', 'FORMULA', 'EXP_MMW', 'SCORE', 'OC']])
```

## Common Scenarios

### Scenario 1: High-Resolution Data with Tight Mass Accuracy

```python
annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    formula_flp=formula,
    transformation_list=trans_list,
    ntheo=12,          # More peaks
    min_overlap=0.7,   # Stricter
    max_msigma=3,      # Lower threshold
    max_mmw_ppm=5,     # Tight tolerance
    baseline=1000
)
```

### Scenario 2: Noisy Data or Low Resolution

```python
annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    formula_flp=formula,
    transformation_list=trans_list,
    ntheo=6,           # Fewer peaks
    min_overlap=0.4,   # More permissive
    max_msigma=15,     # Higher threshold
    max_mmw_ppm=20,    # Loose tolerance
    baseline=100
)
```

### Scenario 3: MS/MS Fragment Analysis

```python
# Process with MSMS=True
result = process_scan(
    test_scan=raw_data,
    MSMS=True,
    baseline=50,
    min_mw=0,
    max_mw=2000
)

# Use molecular database instead of transformations
fragment_db = pd.DataFrame({
    'CPD': ['y12', 'y11', 'w8'],
    'FORMULA': ['C10H13N2O8P', 'C9H12N2O7P', 'C8H11N2O6P']
})

annotated = annotate_scan_targeted(
    scan_processed_aggregated=result['scan_processed_aggregated'],
    formula_flp="",
    cpd_flp="",
    mdb=fragment_db,
    ntheo=6,
    min_overlap=0.4,
    baseline=50
)
```

### Scenario 4: Large Transformation List

```python
# For large lists, consider filtering first
large_trans_list = pd.read_csv('huge_transformation_list.txt', sep='\t')

# Filter by mass range
filtered = large_trans_list[
    (large_trans_list['Delta.AVG.MW'].abs() < 500)
]

annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    formula_flp=formula,
    transformation_list=filtered,
    ...
)
```

## Key Functions and What They Do

| Function | Purpose | Example |
|----------|---------|---------|
| `parse_formula()` | Parse formula string | `parse_formula("C12H22O11")` → `{'C': 12, 'H': 22, 'O': 11}` |
| `calculate_mass()` | Calculate mono/avg mass | `calculate_mass(formula, 'mono')` |
| `get_isotope_distribution()` | Get isotope pattern | `masses, ints = get_isotope_distribution(formula, ntheo=10)` |
| `cut_mmw_list()` | Cluster peaks | `clusters = cut_mmw_list(mw_array, int_array, window=10)` |
| `expand_transformation_list()` | Add formulas to trans list | `ifl, trans = expand_transformation_list(formula, trans_list)` |
| `annotate_scan_targeted()` | Main annotation function | (See examples above) |

## Output Columns Explained

### Annotated Spectrum (`scan` table)

- **MW**: Molecular weight of peak
- **Response**: Intensity
- **CPD**: Matched compound name
- **FORMULA**: Matched elemental formula
- **SCORE**: Chi-square score (lower = better match)
- **OC_SCORE**: Overlap coefficient (0-1, higher = better)
- **COEF**: Mixture coefficients (coef1:coef2)

### Feature Summary (`feature` table)

- **FEATURE**: Unique feature identifier
- **FORMULA**: Elemental formula
- **THEO_MMW**: Theoretical monoisotopic mass
- **THEO_AMW**: Theoretical average mass
- **EXP_MMW**: Experimental monoisotopic mass
- **EXP_AMW**: Experimental average mass
- **EXP_MMW_PPM**: PPM error for monoisotopic mass
- **EXP_AMW_DEV**: Absolute error for average mass
- **SCORE**: Chi-square score
- **OC**: Overlap coefficient
- **RESPONSE**: Total intensity of feature

## Troubleshooting

### Problem: No features found

**Possible causes:**
- Transformation list doesn't match sample
- Baseline threshold too high
- Mass accuracy too strict (`max_mmw_ppm` too small)
- Isotope pattern overlap too strict (`min_overlap` too high)

**Solutions:**
```python
# Try relaxing parameters
annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    transformation_list=trans_list,
    max_mmw_ppm=20,   # Increase from 10
    min_overlap=0.4,  # Decrease from 0.6
    max_msigma=10,    # Increase from 5
    baseline=500      # Decrease if your data is noisy
)
```

### Problem: Too many false positives

**Possible causes:**
- Parameters too permissive
- Baseline too low
- Isotope overlap threshold too low

**Solutions:**
```python
# Try tightening parameters
annotated = annotate_scan_targeted(
    scan_processed_aggregated=aggregated,
    transformation_list=trans_list,
    max_mmw_ppm=5,    # Decrease from 10
    min_overlap=0.7,  # Increase from 0.6
    max_msigma=3,     # Decrease from 5
    baseline=1000     # Increase
)
```

### Problem: Memory error

**For large datasets:**
```python
# Process in chunks
chunk_size = 1000
for i in range(0, len(data), chunk_size):
    chunk = data.iloc[i:i+chunk_size]
    result = annotate_scan_targeted(chunk, ...)
    # Save to file to free memory
    result.to_csv(f'result_{i}.csv')
```

## Running Tests

```bash
# Run all tests
pytest test_annotate_scan.py -v

# Run specific test class
pytest test_annotate_scan.py::TestFormulaParser -v

# Run with coverage
pytest test_annotate_scan.py --cov=annotate_scan
```

## Performance Tips

1. **Filter transformation list** by mass range before processing
2. **Reduce `ntheo`** for faster calculations if resolution allows
3. **Process large datasets** in chunks
4. **Cache isotope distributions** for repeated formulas
5. **Use lower baseline** for noisy data to avoid too many candidates

## Related Functions in process_scan.py

You'll also need:

```python
from process_scan import process_scan

# Main function to prepare data
result = process_scan(
    test_scan=raw_ms_data,
    polarity='Negative',
    MSMS=False,
    baseline=1000,
    min_charge=3,
    max_charge=12,
    min_mz=500,
    max_mz=1500,
    min_mw=4000,
    max_mw=12000,
    mz_error=0.02
)
```

## Resources

- R Source: https://github.com/daniellyz/OligoDistiller/blob/master/R/annotate_scan_targeted.R
- Example data: https://github.com/daniellyz/MESSAR/tree/master/MESSAR_WEBSERVER_DEMO
- brainpy documentation: https://github.com/mobiusklein/brainpy
- Mass spectrometry basics: https://www.upress.psu.edu/books/9780271033600-proteomics/

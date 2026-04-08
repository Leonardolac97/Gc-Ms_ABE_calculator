# GC-MS ABE Calculator

A Python desktop application for processing GC-MS `.csv` files from ABE mixture reaction experiments. The tool loads a classified peak table, filters product entries marked as `KEEP`, excludes reactants from the selectivity basis, and generates selectivity and conversion metrics.

---

## Expected input file structure

The input CSV must follow a structured peak table format (e.g., exported from GC–MS software).

### Required columns

| Column Name       | Description |
|------------------|------------|
| RT               | Retention time (min) |
| Compound Name    | Identified compound name |
| Formula          | Molecular formula (e.g., C7H14O) |
| Carbon Number    | Optional; inferred from formula if missing |
| Area %           | Relative peak area (%) |
| Classification   | Must include `KEEP`, `DISCARD`, or `FLAG` |
| Notes            | Optional annotation |
| Family           | Chemical family (e.g., Alcohols, Ketones, Reactants) |

### Example input (excerpt)

| RT   | Compound Name      | Formula | Carbon Number | Area % | Classification | Notes                  | Family     |
|------|------------------|---------|---------------|--------|----------------|------------------------|-----------|
| 2.91 | Ethanol           | C2H6O   | 2             | 12.3   | KEEP           | Reactant               | Reactants |
| 3.01 | Acetone           | C3H6O   | 3             | 10.1   | KEEP           | Reactant               | Reactants |
| 3.14 | Isopropanol       | C3H8O   | 3             | 4.2    | KEEP           | Transfer hydrogenation | Alcohols  |
| 5.13 | 1-Butanol         | C4H10O  | 4             | 8.5    | KEEP           | Reactant               | Reactants |

### Important rules

- Only rows with `Classification == KEEP` are used  
- Rows with `Family == Reactants` are **excluded** from selectivity calculations  
- `Area %` is used as the calculation basis  
- Carbon number is inferred automatically if not provided  

---

## Data origin

The input files are assumed to be exported from GC–MS analysis software (e.g., Agilent ChemStation or MassHunter) and manually classified into:

- `KEEP` → valid products  
- `DISCARD` → noise or irrelevant peaks  
- `FLAG` → optional review  

The workflow is designed for catalytic ABE upgrading experiments.

---

## Features

- Load sample CSV files  
- Optional baseline (`0 h`) for conversion calculations  
- Manual A0 input for flexibility  
- Automatic export of processed data  
- Interactive plots using PyQt6 and Matplotlib  

---

## Installation

### Option 1 — Conda

```bash
conda env create -f environment.yml
conda activate gc-ms-abe-calculator
python GC_MS_ABE_calculator.py
```

### Option 2 — pip

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate

pip install -r requirements.txt
python GC_MS_ABE_calculator.py
```

## Output files

After loading a sample, the application creates a results folder next to the input file:

- `keep_products_table.csv`
- `family_selectivity.csv`
- `carbon_selectivity.csv`
- `carbon_bins_selectivity.csv`
- `product_selectivity.csv`
- `conversion.csv` (only when conversion can be calculated)


## Citation

Please cite this software using the metadata provided in `CITATION.cff`.

## License

This project is released under the MIT License. See `LICENSE`.

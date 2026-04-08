# GC-MS ABE Calculator

A Python desktop application for processing GC-MS `.csv` files from ABE mixture reaction experiments. The tool loads a classified peak table, filters product entries marked as `KEEP`, excludes reactants from the selectivity basis, and generates:

- family selectivity
- carbon-number selectivity
- carbon-bin selectivity (`C2–C6` and `C7–C15`)
- selected-product selectivity
- reactant conversion, when a baseline (`0 h`) file or manual initial area percentages are provided

The graphical user interface is built with PyQt6 and the plots are generated with Matplotlib.

## Features

- Loads a sample `.csv` file and processes the data automatically
- Optional baseline (`0 h`) import for conversion calculations
- Optional manual entry of initial reactant area percentages
- Exports processed tables as `.csv`
- Creates plots in separate tabs for:
  - family selectivity
  - carbon number
  - carbon bins
  - selected products
  - conversion
- Opens the results folder directly from the interface

## Expected input columns

The sample CSV file must contain the following columns:

- `RT`
- `Compound Name`
- `Formula`
- `Area %`
- `Classification`
- `Notes`
- `Family`

If `Carbon Number` is not present, it is inferred from the molecular formula.

## How it works

1. The application reads the input CSV.
2. Only rows with `Classification == KEEP` are retained.
3. Rows with `Family == Reactants` are excluded from the product selectivity basis.
4. Selectivity is calculated from the remaining `Area %` values.
5. Reactant conversion is computed from:
   - a baseline file (`0 h`), or
   - manually entered initial area percentages (`A0`)

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

## Notes for publication

Before publishing the repository, it is recommended to:

- remove personal or internal-only comments that are not needed for public release
- confirm that the target product aliases match the intended chemistry scope
- add a short example input file or screenshot if redistribution is allowed
- state the license clearly in the repository root

## Citation

Please cite this software using the metadata provided in `CITATION.cff`.

## License

This project is released under the MIT License. See `LICENSE`.

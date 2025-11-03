# Python Simulation Scripts

This folder contains the Python scripts for cancer progression simulations.

## Files

- `cancer_sim.py` - Cancer progression simulator with therapy and parameter variability
- `cancer_sim2.py` - Stage I Cancer Model with multi-therapy schedules (surgery, chemo, radiotherapy)
- `Lung_cancer.py` - Lung cancer specific simulation
- `import json.py` - JSON data import/export utilities
- `requirements.txt` - Python dependencies

## Setup

1. Create and activate a virtual environment:

```bash
cd python
python3 -m venv venv
source venv/bin/activate  # On macOS/Linux
# or
venv\Scripts\activate  # On Windows
```

2. Install dependencies:

```bash
pip install -r requirements.txt
```

## Running Simulations

### Run cancer_sim.py:

```bash
python cancer_sim.py
```

This generates a plot showing cancer progression with chemo and radiotherapy.

### Run cancer_sim2.py:

```bash
python cancer_sim2.py
```

This generates multiple output files in `../outputs/`:

- `pulse_chemo.png`
- `metronomic_chemo.png`
- `weekly_chemo.png`
- `radiotherapy.png`
- `stage1_surgery_repop_outputs.json`

The JSON output file is used by the Next.js frontend application.

## Output

Simulation outputs are saved to `../outputs/` directory (one level up from this folder).

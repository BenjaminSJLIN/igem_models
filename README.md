# iGEM Model Pipeline

This repository contains the mathematical models and simulation pipeline for the iGEM project.

## Model 0.5: CRISPR-Cas12a Activation Dynamics

This module simulates the expression dynamics of the Cas12a protein activated by L-Arabinose (Ara). It utilizes a system of Ordinary Differential Equations (ODEs) to model both mRNA and Protein levels over time based on the Hill Equation.

### Project Structure

- `data/`: Contains raw, intermediate, and final data outputs.
- `models/`: Contains the logic and equations for each stage.
- `outputs/`: Generated plots and visualizations.
- `utils/`: Shared utilities for plotting and data processing.

## Setup Instructions

1. **Create a virtual environment:**
   ```bash
   python -m venv venv
   ```

2. **Activate the virtual environment:**
   - **Windows (PowerShell):** `.\venv\Scripts\Activate.ps1`
   - **Windows (CMD):** `.\venv\Scripts\activate.bat`
   - **Mac/Linux:** `source venv/bin/activate`

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Running the Simulation

Execute the main orchestrator script:
```bash
python main_pipeline.py
```

This will:
1. Load parameters from `config.yaml`.
2. Run the ODE simulation across different Arabinose concentrations.
3. Save the simulated data to `data/02_intermediate/model_0_5_result.csv`.
4. Generate and save time-series and dose-response plots in `outputs/model_0_5_plots/`.

## Configuration

You can easily adjust biological parameters (e.g., maximum transcription rate, degradation rates, Hill coefficient) in the `config.yaml` file without modifying the code.

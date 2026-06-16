# MT2D DirectJ-Broyden Demo

This repository provides a MATLAB demo for 2-D magnetotelluric joint TE-TM inversion using fixed synthetic data and cached direct-sensitivity Jacobians.

## Requirements

- MATLAB R2023a or newer is recommended.
- No additional MATLAB toolbox is required for the provided demo workflow.

## Repository structure

```text
MT2D-DirectJ-Broyden-demo/
├── data/
│   ├── ascii_layered_100_10/
│   ├── ascii_vertical_100_10/
│   └── jacobian/
├── outputs/
├── src/
├── CITATION.cff
├── EXPECTED_RESULTS.txt
├── LICENSE
├── README.md
└── run_full_inversion_layered_100_10.m
```

## How to run the inversion demo

1. Open MATLAB.
2. Set the repository root as the current folder:

```matlab
cd('path_to/MT2D-DirectJ-Broyden-demo')
```

3. Run:

```matlab
run_full_inversion_layered_100_10
```

The script loads the selected synthetic data folder, loads the cached direct-sensitivity Jacobians, runs the joint TE-TM inversion with controlled Broyden updates, prints the RMS and timing summary in the MATLAB Command Window, and displays the model, RMS convergence, and data-fit figures.

## Selecting a synthetic model

The synthetic case is selected in:

```text
src/Inversi_TETM_DirectJ_Broyden.m
```

Edit the line:

```matlab
model_type = 'layered_100_10';
```

Available cases in this repository:

```matlab
model_type = 'layered_100_10';
model_type = 'vertical_100_10';
```

Each case requires a matching input folder under `data/` and matching cached Jacobian files under `data/jacobian/`.

For example, `model_type = 'layered_100_10'` uses:

```text
data/ascii_layered_100_10/
data/jacobian/J_TE_direct_shared_layered_100_10.mat
data/jacobian/J_TM_direct_shared_layered_100_10.mat
```

and `model_type = 'vertical_100_10'` uses:

```text
data/ascii_vertical_100_10/
data/jacobian/J_TE_direct_shared_vertical_100_10.mat
data/jacobian/J_TM_direct_shared_vertical_100_10.mat
```

## Main settings

The release demo uses cached Jacobians by default:

```matlab
force_build_initial_jacobian = false;
load_only_jacobian           = true;
```

Figure and output options are controlled by:

```matlab
show_figures = true;
save_results = false;
save_figures = false;
```

Set `save_results = true` only if you want to save a compact `.mat` result file.

Set `save_figures = true` only if you want MATLAB to save the displayed figures.

## Expected output

For the `layered_100_10` case, the expected result is approximately:

```text
Initial joint RMS : 7.44
Final joint RMS   : 1.08
Final TE RMS      : 1.07
Final TM RMS      : 1.09
```

Small numerical differences may occur because of MATLAB version, sparse linear solver behaviour, operating system, and hardware-dependent numerical precision.

## Citation

If you use this code, please cite the associated manuscript and repository. Citation metadata are provided in `CITATION.cff`.

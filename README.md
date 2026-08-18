<img src="Hyperelastic_Materials.png" alt="Description" style="max-width: 100%; height: auto;">


# Installation requirements

To install the required dependencies, run the following command in your terminal:

```bash
pip install matplotlib scipy sympy numpy
```


# 🏛️ The Foundations of Nonlinear Solid Mechanics

This repository presents a compact introduction to constitutive modeling for hyperelastic materials such as rubbers, soft polymers, and biological tissues. These materials can undergo large deformations, so their mechanical response must be described with nonlinear models rather than classical linear assumptions.

The code in this repository is organized around three main ideas:

- describing the kinematics of deformation,
- developing constitutive models,
- enforcing equilibrium conditions.

It provides a practical toolkit for studying hyperelastic behavior through benchmark deformation cases such as uniaxial tension and pure shear. The scripts use Python and `SymPy` to carry out symbolic calculations in a clear and reproducible way.

This repository is intended as an introductory guide for anyone who wants to explore constitutive modeling for nonlinear solids. It includes implementations of several well-known hyperelastic models, including:

- **Neo-Hookean**:


$$
W_{\mathrm{NH}} = C_{10}(I_1-3)
$$


- **Mooney–Rivlin**:


$$
W_{\mathrm{MR}} = C_{10}(I_1 - 3) + C_{01}(I_2 - 3)
$$


- **Gent**:

$$
W_{\mathrm{Gent}} = -\frac{\mu J_m}{2}\ln\!\left(1 - \frac{I_1 - 3}{J_m}\right)
$$

- **Yeoh**:

$$
W_{\mathrm{Yeoh}} = C_{10}(I_1-3) + C_{20}(I_1-3)^2 + C_{30}(I_1-3)^3.
$$

- **Anssari-Benam**:

$$

W({I}_1,{I}_2) = 3^{\frac{n_1-1}{2n_1}} \mu_1 N_1 \left[\frac{1}{3N_1(n_1-1)} ({I}_1-3)^{\beta_1} - \ln\left[ \left( \frac{{I}_1-3N_1}{3-3N_1} \right)^{\beta_1} \right] + c_{20} \left[ \left(\frac{{I}_2}{3}\right)^{\varepsilon_1} -1 \right] \right]

$$

These models are commonly used to represent the behavior of soft materials under large deformation. Before they can be used in analysis, their material parameters must be fitted to experimental data for the material of interest.

## What this repository provides

- A clear introduction to hyperelastic constitutive modeling
- Symbolic derivations for common benchmark problems
- Python-based scripts for model evaluation
- A starting point for fitting material constants to test data
- A foundation for future finite element or research applications

## Typical use cases

- studying the response of soft materials under load
- comparing different hyperelastic material models
- generating stress–stretch curves for benchmark tests
- supporting research and teaching in nonlinear mechanics

<table>
  <tr>
    <td width="100%">
      <img src="Hyperelastic_Materials2.jpg" alt="Figure 2" style="border-radius: 8px; width: 100%;">
    </td>
  </tr>
</table>

---


# 📉 Calibration
The `calibration` folder contains a standard raw Python code as well as a Python app (more user-friendly) for fitting well-known material models (mentioned before).
The code uses experimental data for uniaxial tension and pure shear tests.

## Calibration raw code

In `Calibration_raw.py`, experimental data should be provided as a .txt file named `test_data.txt`. The code jointly calibrates the chosen models to the experimental data using nonlinear least-squares optimization. The first Piola–Kirchhoff stress is used in this case.
At the end, the material parameters are estimated and reported. The quality of calibration is evaluated using metrics:

- Relative root-mean-square error
- Root-mean-square error (**RMSE**)
- Normalized root-mean-square error (**NRMSE**)
- Coefficient of determination ($R^2$)
- Maximum absolute error

The fitted material constants and corresponding goodness-of-fit measures are then reported for each constitutive model.
Finally, the code generates comparison plots of the calibrated model predictions against the experimental data, saves the figures in `.svg` and `.tiff` formats, and writes a formatted summary of the calibration results to `fit_results.txt`.

> **Important:** Appropriate parameter bounds and initial guesses are critical for obtaining stable and physically meaningful solutions from nonlinear least-squares optimization.


## Calibration App

A Python app for calibrating hyperelastic constitutive models to experimental stress-stretch data saved as `Calibration_App`.

## What it does

Loads uniaxial tension and pure shear test data from a CSV file, fits one of five hyperelastic models to the data using nonlinear least squares, and displays the fitted parameters alongside comparison plots.

## Supported Models

| Model | Parameters |
|---|---|
| Neo-Hookean | C₁₀ |
| Mooney-Rivlin | C₁₀, C₀₁ |
| Gent | μ, Jₘ |
| Yeoh | C₁, C₂, C₃ |
| Anssari-Benam | μ, N, n, β, C₂₀, ε |

## Input Format

A `.csv` file with four columns (no header required):


---


<table>
  <tr>
    <td width="100%">
        <img src="calibration_fit.jpg" alt="Description" style="max-width: 100%; height: auto;">
    </td>
  </tr>
</table>

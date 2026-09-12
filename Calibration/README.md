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






### Data Calibration Utility (`Calibration_App.py`)

The `Calibration_App.py` script is designed to calibrate hyperelastic constitutive models against experimental data.
It employs a nonlinear least-squares optimization algorithm to jointly calibrate the models against both uniaxial and pure shear test results.

#### Input Requirements
The application requires a `.csv` file containing experimental data from uniaxial and pure shear tests. The data must be structured with the following column headers:

* `stretch_uniaxial`
* `stress_uniaxial`
* `stretch_pure_shear`
* `stress_pure_shear`

#### Implemented Constitutive Models
The code includes the following hyperelastic material models:

* Neo-Hookean
* Mooney-Rivlin
* Yeoh
* Gent
* Anssari-Benam


| Model | Parameters |
|:---:|:---:|
| Neo-Hookean | C₁₀ |
| Mooney-Rivlin | C₁₀, C₀₁ |
| Gent | μ, Jₘ |
| Yeoh | C₁, C₂, C₃ |
| Anssari-Benam | μ, N, n, β, C₂₀, ε |


#### Output and Visualization
Upon execution, the application provides the following outputs:

1.  **Visualization:** Calibration curves, goodness-of-fit metrics (e.g., $R^2$, RMSE), and optimized material parameters are rendered in dedicated UI panels for immediate review.
2.  **Exported Files:** 
    *   A vector graphic of the calibration plots (`.svg`).
    *   A summary report containing fitting metrics and optimized material parameters (`.txt`), saved automatically to the working directory.



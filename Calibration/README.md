### Data Calibration Utility (`Calibration_App.py`)

The `Calibration_App.py` script is designed to calibrate hyperelastic constitutive models against experimental data.

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



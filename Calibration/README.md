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





This is a concise, professional `README.md` tailored for `grabit.py`, designed to get a researcher up and running quickly.

***

# Image Graph Digitizer

**Image Graph Digitizer** is a lightweight, Python-based GUI tool designed for extracting numerical data from image-based plots. It provides a simple calibration workflow and allows for flexible data export, making it a useful utility for researchers needing to digitize data from legacy plots, paper figures, or experimental images.


<table>
  <tr>
    <td width="100%">
      <img src="Calibration\Digitizer_GUI.png" alt="Figure 1" style="border-radius: 8px; width: 100%;">
    </td>
  </tr>
</table>

## Features
*   **Intuitive Calibration:** Supports linear and logarithmic axes (X and Y).
*   **Navigation:** Zoom and Pan capabilities for precise point selection.
*   **Dataset Management:** Create, rename, and organize multiple datasets within a single session.
*   **Flexible Exports:** Save data in various formats including `.txt`, `.csv`, `.npz` (NumPy), and `.mat` (MATLAB).

## Requirements
*   **Python 3.x**
*   **Pillow** (Required for image processing): `pip install Pillow`
*   **NumPy & SciPy** (Optional, for advanced export formats): `pip install numpy scipy`

## Quick Start
1.  **Launch:** Run the script using `python grabit.py`.
2.  **Load:** Click "Load Image..." to import your plot.
3.  **Calibrate (Crucial):** 
    *   Click "Calibrate".
    *   Follow the status bar prompts to click four points on the graph: **X-Origin**, **X-Max**, **Y-Origin**, and **Y-Max**.
    *   Enter the corresponding real-world values for these points when prompted.
4.  **Grab:** Click "Grab Points" and begin clicking on your data series.
5.  **Finish:** Press **Enter** when done to finalize the dataset, then save/export it using the sidebar.

## Keyboard Shortcuts
*   **`a` / `b`**: Zoom In / Zoom Out.
*   **`Space`**: Fit image to window.
*   **`r`**: Restore zoom level.
*   **`Backspace` / `Delete`**: Remove the last grabbed point (during "Grabbing" mode).
*   **`Enter`**: Complete point acquisition.
*   **Mouse Middle-Click (Drag)**: Pan the image.


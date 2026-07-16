<img src="Hyperelastic_Materials.png" alt="Description" style="max-width: 100%; height: auto;">

🏛️ The Foundations of Nonlinear Solid Mechanics
Constitutive modeling is the essential mathematical framework that enables us to predict how hyperelastic materials—such as rubbers, bio-tissues, and soft polymers—respond to external loads by defining a Strain Energy Density Function ($W$) 📐. Unlike linear materials, hyperelastic solids exhibit extreme geometric and material nonlinearities, often undergoing massive deformations 🌀. By accurately capturing this complex mechanical behavior, constitutive models allow engineers to perform reliable Finite Element Analysis (FEA), which is indispensable for modern engineering design 🏎️.



To put these theories into practice, the GitHub repository below serves as a comprehensive toolkit, offering derived equations for fundamental benchmark problems such as uniaxial tension and pure shear 🛠️. Within the codebase, you will quickly discover the power and convenience of the SymPy library as a versatile tool for symbolic calculations in Python 🐍. These scripts are designed to act as an introductory tutorial, providing a clear pathway to dive into the world of constitutive modeling for nonlinear solids. High-fidelity implementations of well-known models, including the Neo-Hookean, Mooney-Rivlin, and Yeoh formulations, are all developed within the code to jumpstart your research 🌐.

This repository provides comprehensive information regarding **hyperelastic materials** and their **constitutive modeling**.

Based on the available data, equations have been derived for **uniaxial tension** and **pure shear** problems. 

In the `extension-torsion` file, the components of the **Cauchy stress tensor** are provided for the problem of **extension superimposed on torsion**, utilizing a well-known strain energy function. 

---

### 📉 Calibration

The `calibration` folder contains standard code for calibrating hyperelastic materials using several well-known constitutive models, including:
* **Mooney-Rivlin**
* **Yeoh**
The code experimental data from test_data.txt, then jointly calibrates three hyperelastic constitutive models—Mooney-Rivlin, Yeoh, and a custom Anssari-Benam-type model—by fitting their nominal stress predictions to both datasets using nonlinear least-squares optimization. It defines the stress response equations (First Pila P) for each loading mode, estimates the material parameters, and computes error metrics such as relative RMSE, RMSE, NRMSE, R2, and maximum error, and reports the fitted constants and goodness-of-fit for each model. 
> prescribed bounds and initial points are needed (critical) for least-square optimization.
Finally, it generates comparison plots of the calibrated model curves against the experimental data, saves the figures as .svg and .tiff, and writes a formatted summary of the calibration results to fit_results.txt.


---

### 🐍 Dependencies & Libraries

Regarding the Python files, the following libraries have been utilized:
* `numpy`
* `pandas`
* `matplotlib`
*  `sympy` 
*  `scipy`


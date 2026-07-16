"""
Code by Mohammad Ali Safaei

mohsafaei.github.io
email: mohammadsf1998@gmail.com

Experimental data source: https://doi.org/10.1007/s10659-024-10055-y 
this code nees to be in the same folder with test_data.txt as the input file.
The Mooney-Rivlin, Yeoh, and Anssari-Benam models are simultaneously calibrated using uniaxial tension and pure shear data. 
Plots compare the calibrated constitutive models with experimental data points.
"""

# -----------------------------------------------------------------------
# 0. Import libraries
# -----------------------------------------------------------------------

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------
# 1. Load experimental data from the txt file
# -----------------------------------------------------------------------

data_path = "test_data.txt"

namespace = {"np": np}
with open(data_path, "r") as f:
    exec(f.read(), namespace)

uniaxial_data = namespace["uniaxial_data"]
shear_data = namespace["shear_data"]

lam_u, P_u = uniaxial_data[1:, 0], uniaxial_data[1:, 1]   # drop trivial (1,0) point
lam_s, P_s = shear_data[1:, 0], shear_data[1:, 1]

#print(f"Uniaxial tension data:  n={len(lam_u)}, stretch range [{lam_u.min():.2f}, {lam_u.max():.2f}]")
#print(f"Pure shear data:        n={len(lam_s)}, stretch range [{lam_s.min():.2f}, {lam_s.max():.2f}]")

# -----------------------------------------------------------------------
# 2. Model stress functions (nominal / first Piola-Kirchhoff stress)
# -----------------------------------------------------------------------

def mr_uniaxial(lam, c10, c01):
    return 2*c10*(lam - 1/lam**2) + 2*c01*(1 - 1/lam**3)   #correct

def mr_pureshear(lam, c10, c01):
    # I1 = I2 for pure shear -> stress simplifies to 2*(lam-1/lam^3)*(c10+c01)
    return 2*(lam - 1/lam**3) * (c10 + c01)

def I1_uniaxial(lam):
    return lam**2 + 2/lam

def I1_pureshear(lam):
    return lam**2 + 1/lam**2 + 1


def yeoh_uniaxial(lam, c1, c2, c3):
    I1 = I1_uniaxial(lam)
    dW = c1 + 2*c2*(I1 - 3) + 3*c3*(I1 - 3)**2
    return 2*(lam - 1/lam**2) * dW

def yeoh_pureshear(lam, c1, c2, c3):
    I1 = I1_pureshear(lam)
    dW = c1 + 2*c2*(I1 - 3) + 3*c3*(I1 - 3)**2
    return 2*(lam - 1/lam**3) * dW

#-------------------------------------------
def dW_dI1(I1, mu1, N1, n1, beta1):
    A = 3*(n1-1)*mu1*N1 / (2*n1)
    return A * beta1 * ((I1-3)**(beta1-1) / (3*N1*(n1-1)) - 1/(I1-3*N1))

def dW_dI2(I2, C20, eps1):
    return C20 * eps1 / 3 * (I2/3)**(eps1-1)

def custom_uniaxial(lam, mu1, N1, n1, beta1, C20, eps1):
    I1 = lam**2 + 2/lam
    I2 = 2*lam + 1/lam**2
    w1 = dW_dI1(I1, mu1, N1, n1, beta1)
    w2 = dW_dI2(I2, C20, eps1)
    return 2*(lam - 1/lam**2)*w1 + 2*(1 - 1/lam**3)*w2

def custom_pureshear(lam, mu1, N1, n1, beta1, C20, eps1):
    I1 = lam**2 + 1 + 1/lam**2
    I2 = I1  # I1 == I2 in pure shear
    w1 = dW_dI1(I1, mu1, N1, n1, beta1)
    w2 = dW_dI2(I2, C20, eps1)
    return 2*(lam - 1/lam**3)*(w1 + w2)
# 5. Joint calibration: Custom model
def custom_residuals_joint(x):
    mu1, N1, n1, beta1, C20, eps1 = x
    r_u = P_u - custom_uniaxial(lam_u, mu1, N1, n1, beta1, C20, eps1)
    r_s = P_s - custom_pureshear(lam_s, mu1, N1, n1, beta1, C20, eps1)
    return np.concatenate([r_u, r_s])

fit_custom = least_squares(
    custom_residuals_joint,
    x0=[0.05, 0.9, 0.7, 2.5, 2.5, 0.25],
    bounds=(
        [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6],  # lower
        [10.0, 1.0, 10.0,  5.0, 10.0, 5.0 ]   # upper
    )
)
mu1, N1, n1, beta1, C20, eps1 = fit_custom.x

#----------------------------------------

def rel_rmse(y_exp, y_model):
    rmse = np.sqrt(np.mean((y_exp - y_model)**2))
    return rmse / np.mean(np.abs(y_exp)) * 100

# -----------------------------------------------------------------------
# 3. Joint calibration: Mooney-Rivlin & Yeoh
# -----------------------------------------------------------------------
def mr_residuals_joint(x):
    c10, c01 = x
    P_u_8 = P_u[:8]
    lam_u_8 = lam_u[:8]
    r_u = P_u_8 - mr_uniaxial(lam_u_8, c10, c01)
    r_s = P_s - mr_pureshear(lam_s, c10, c01)
    return np.concatenate([r_u, r_s])

fit_mr = least_squares(
    mr_residuals_joint,
    x0=[0.1, 0.05],          # avoid starting exactly equal too
    bounds=([1e-6, 1e-3], [5, 5])   # c01 lower bound raised to 1e-3
)
c10, c01 = fit_mr.x


def yeoh_residuals_joint(x):
    c1, c2, c3 = x
    r_u = P_u - yeoh_uniaxial(lam_u, c1, c2, c3)
    r_s = P_s - yeoh_pureshear(lam_s, c1, c2, c3)
    return np.concatenate([r_u, r_s])

fit_yeoh = least_squares(yeoh_residuals_joint, x0=[0.2, -0.015, 0.003],
                          bounds=([-1, -1, -1], [1, 1, 1]))
c1, c2, c3 = fit_yeoh.x

# -----------------------------------------------------------------------
# 5. Report calibrated constants and errors
# -----------------------------------------------------------------------

print("\n=== Mooney-Rivlin: joint calibration (uniaxial + pure shear) ===")
print(f"c10 = {c10:.4f} MPa,  c01 = {c01:.4f} MPa")
print(f"rel. RMSE (uniaxial)   = {rel_rmse(P_u, mr_uniaxial(lam_u, c10, c01)):.2f} %")
print(f"rel. RMSE (pure shear) = {rel_rmse(P_s, mr_pureshear(lam_s, c10, c01)):.2f} %")

print("\n=== Yeoh: joint calibration (uniaxial + pure shear) ===")
print(f"c1 = {c1:.5f} MPa,  c2 = {c2:.5f} MPa,  c3 = {c3:.5f} MPa")
print(f"rel. RMSE (uniaxial)   = {rel_rmse(P_u, yeoh_uniaxial(lam_u, c1, c2, c3)):.2f} %")
print(f"rel. RMSE (pure shear) = {rel_rmse(P_s, yeoh_pureshear(lam_s, c1, c2, c3)):.2f} %")


print("\n=== Custom model: joint calibration (uniaxial + pure shear) ===")
print(f"mu1={mu1:.4f}, N1={N1:.4f}, n1={n1:.4f}, beta1={beta1:.4f}, C20={C20:.4f}, eps1={eps1:.4f}")
print(f"rel. RMSE (uniaxial)   = {rel_rmse(P_u, custom_uniaxial(lam_u, mu1, N1, n1, beta1, C20, eps1)):.2f} %")
print(f"rel. RMSE (pure shear) = {rel_rmse(P_s, custom_pureshear(lam_s, mu1, N1, n1, beta1, C20, eps1)):.2f} %")


# -----------------------------------------------------------------------
# 6. Plot calibrated models against the experimental data points
# -----------------------------------------------------------------------

lam_u_fine = np.linspace(lam_u.min(), lam_u.max(), 200)
lam_s_fine = np.linspace(lam_s.min(), lam_s.max(), 200)

fig, axes = plt.subplots(3, 2, figsize=(9., 12.0), constrained_layout=True)
(ax_mr_u, ax_mr_s), (ax_yh_u, ax_yh_s),(ax_ab_u, ax_ab_s) = axes

# --- Mooney-Rivlin: Uniaxial tension ---
ax_mr_u.scatter(lam_u, P_u, marker='h',facecolor='gray',edgecolor='black', label='Uniaxial tension data', zorder=3)
ax_mr_u.plot(lam_u_fine, mr_uniaxial(lam_u_fine, c10, c01),
             linestyle='-', color='darkorange', linewidth=3.5, label='Mooney-Rivlin (joint fit)')
ax_mr_u.set_xlabel(r'$\lambda$', fontsize=13)
ax_mr_u.set_ylabel(r'P [MPa]', fontsize=13)
ax_mr_u.legend()
ax_mr_u.set_xlim(1.0, 3)
ax_mr_u.set_ylim(0.0, 2.5)
ax_mr_u.grid(True, which='major', linestyle='--', alpha=0.4)
ax_mr_u.tick_params(axis='both', labelsize=11, direction='in', pad=8)
ax_mr_u.autoscale(False)

# --- Mooney-Rivlin: Pure shear ---
ax_mr_s.scatter(lam_s, P_s, marker='h',facecolor='gray',edgecolor='black', label='Pure shear data', zorder=3)
ax_mr_s.plot(lam_s_fine, mr_pureshear(lam_s_fine, c10, c01),
             linestyle='-', color='darkorange', linewidth=3.5, label='Mooney-Rivlin (joint fit)')
ax_mr_s.set_xlabel(r'$\lambda$', fontsize=13)
ax_mr_s.set_ylabel(r'P [MPa]', fontsize=13)
ax_mr_s.legend()
ax_mr_s.set_xlim([1.0, 1.5])
ax_mr_s.set_ylim([0.0, 1.0])
ax_mr_s.grid(True, which='major', linestyle='--', alpha=0.4)
ax_mr_s.tick_params(axis='both', labelsize=11, direction='in', pad=8)

# --- Yeoh: Uniaxial tension ---
ax_yh_u.scatter(lam_u, P_u, marker='h',facecolor='gray',edgecolor='black', label='Uniaxial tension data', zorder=3)
ax_yh_u.plot(lam_u_fine, yeoh_uniaxial(lam_u_fine, c1, c2, c3),
             linestyle='-', color='#ab0505', linewidth=3.5, label='Yeoh (joint fit)')
ax_yh_u.set_xlabel(r'$\lambda$', fontsize=13)
ax_yh_u.set_ylabel(r'P [MPa]', fontsize=13)
ax_yh_u.legend()
ax_yh_u.set_xlim([1.0, 3.0])
ax_yh_u.set_ylim([0.0, 2.5])
ax_yh_u.grid(True, which='major', linestyle='--', alpha=0.4)
ax_yh_u.tick_params(axis='both', labelsize=11, direction='in', pad=8)

# --- Yeoh: Pure shear ---
ax_yh_s.scatter(lam_s, P_s, marker='h',facecolor='gray',edgecolor='black', label='Pure shear data', zorder=3)
ax_yh_s.plot(lam_s_fine, yeoh_pureshear(lam_s_fine, c1, c2, c3),
             linestyle='-', color='#ab0505', linewidth=3.5, label='Yeoh (joint fit)')
ax_yh_s.set_xlim([1.0, 1.5])
ax_yh_s.set_ylim([0.0, 1.0])
ax_yh_s.set_xlabel(r'$\lambda$', fontsize=13)
ax_yh_s.set_ylabel(r'P [MPa]', fontsize=13)
ax_yh_s.legend()
ax_yh_s.grid(True, which='major', linestyle='--', alpha=0.4)
ax_yh_s.tick_params(axis='both', labelsize=11, direction='in', pad=8)

# --- Anssari-Benam: Uniaxial tension ---
ax_ab_u.scatter(lam_u, P_u, marker='h', facecolor='gray',edgecolor='black', label='Uniaxial tension data', zorder=3)
ax_ab_u.plot(lam_u_fine, custom_uniaxial(lam_u_fine, mu1, N1, n1, beta1, C20, eps1),
             linestyle='-', color='#0f8745', linewidth=3.5, label='Anssari-Benam (joint fit)')
ax_ab_u.set_xlim([1.0, 3.0])
ax_ab_u.set_ylim([0.0, 2.5])
ax_ab_u.set_xlabel(r'$\lambda$', fontsize=13)
ax_ab_u.set_ylabel(r'P [MPa]', fontsize=13)
ax_ab_u.legend()
ax_ab_u.grid(True, which='major', linestyle='--', alpha=0.4)
ax_ab_u.tick_params(axis='both', labelsize=11, direction='in', pad=8)

# --- Anssari-Benam: Pure shear ---
ax_ab_s.scatter(lam_s, P_s, marker='h', facecolor='gray',edgecolor='black', label='Pure shear data', zorder=3)
ax_ab_s.plot(lam_s_fine, custom_pureshear(lam_s_fine, mu1, N1, n1, beta1, C20, eps1),
             linestyle='-', color='#0f8745', linewidth=3.5, label='Anssari-Benam (joint fit)')
ax_ab_s.set_xlim([1.0, 1.5])
ax_ab_s.set_ylim([0.0, 1.0])
ax_ab_s.set_xlabel(r'$\lambda$', fontsize=13)
ax_ab_s.set_ylabel(r'P [MPa]', fontsize=13)
ax_ab_s.legend()
ax_ab_s.grid(True, which='major', linestyle='--', alpha=0.4)
ax_ab_s.tick_params(axis='both', labelsize=11, direction='in', pad=8)
fig.savefig("calibration_fit_separated.svg")
fig.savefig("calibration_fit_separated.tiff")

plt.show()

# -----------------------------------------------------------------------
# 5. Report calibrated constants and errors
# -----------------------------------------------------------------------

def metrics(y_exp, y_model):
    res = y_exp - y_model
    rmse = np.sqrt(np.mean(res**2))
    nrmse = rmse / (y_exp.max() - y_exp.min()) * 100
    ss_res = np.sum(res**2)
    ss_tot = np.sum((y_exp - y_exp.mean())**2)
    r2 = 1 - ss_res / ss_tot
    max_err = np.max(np.abs(res))
    return rmse, nrmse, r2, max_err

def build_rows(u_pred, s_pred):
    P_comb, pred_comb = np.concatenate([P_u, P_s]), np.concatenate([u_pred, s_pred])
    rows = []
    for label, y_e, y_m in [("Uniaxial", P_u, u_pred),
                             ("Pure shear", P_s, s_pred),
                             ("Combined", P_comb, pred_comb)]:
        rows.append((label, *metrics(y_e, y_m)))
    return rows

def format_table(title, params_line, rows):
    header = f"{'Test mode':<12}{'RMSE':>12}{'NRMSE %':>12}{'R2':>10}{'Max error':>14}"
    sep = "-" * len(header)
    lines = [title, params_line, sep, header, sep]
    for label, rmse, nrmse, r2, max_err in rows:
        lines.append(f"{label:<12}{rmse:>12.4f}{nrmse:>12.2f}{r2:>10.4f}{max_err:>14.4f}")
    lines.append(sep)
    return "\n".join(lines)

# ------------------------------------------------------

P_u_mr   = mr_uniaxial(lam_u, c10, c01)
P_s_mr   = mr_pureshear(lam_s, c10, c01)
P_u_yeoh = yeoh_uniaxial(lam_u, c1, c2, c3)
P_s_yeoh = yeoh_pureshear(lam_s, c1, c2, c3)
P_u_ab = custom_uniaxial(lam_u, mu1, N1, n1, beta1, C20, eps1)
P_s_ab = custom_pureshear(lam_s, mu1, N1, n1, beta1, C20, eps1)

mr_rows   = build_rows(P_u_mr, P_s_mr)
yeoh_rows = build_rows(P_u_yeoh, P_s_yeoh)
ab_rows = build_rows(P_u_ab, P_s_ab)

mr_table = format_table(
    "=== Mooney-Rivlin: joint calibration (uniaxial + pure shear) ===",
    f"c10 = {c10:.4f} MPa,  c01 = {c01:.4f} MPa",
    mr_rows
)
yeoh_table = format_table(
    "=== Yeoh: joint calibration (uniaxial + pure shear) ===",
    f"c1 = {c1:.5f} MPa,  c2 = {c2:.5f} MPa,  c3 = {c3:.5f} MPa",
    yeoh_rows
)
ab_table = format_table(
    "=== Anssari-Benam: joint calibration (uniaxial + pure shear) ===",
    f"mu1 = {mu1:.5f} MPa, C20= {C20:.5f} MPa, n1 = {n1:.5f}, "
    f"N1 = {N1:.5f},  beta1= {beta1:.5f}, epsilon_1= {eps1:.5f}",
    ab_rows
)

print("\n" + mr_table)
print("\n" + yeoh_table)
print("\n" + ab_table)

results_path = "fit_results.txt"
with open(results_path, "w") as f:
    f.write("Hyperelastic material calibration - fit results\n")
    f.write("=" * 60 + "\n\n")
    f.write(mr_table + "\n\n")
    f.write(yeoh_table + "\n")
    f.write(ab_table + "\n")

print(f"\nSaved results table: {results_path}")


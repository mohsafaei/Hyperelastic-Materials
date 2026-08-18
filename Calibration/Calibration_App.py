"""
Code by Mohammad Ali Safaei

website:
       mohsafaei.github.io
email: 
       mohammadsf1998@gmail.com
       ma.safaei@ut.ac.ir
Please cite: Multiaxial Finite Strain Behavior of Polydomain Liquid Crystal Elastomers

       
Experimental data source: https://doi.org/10.1007/s10659-024-10055-y 
this code needs to be in the same folder with test_data.txt as the input file.
The Neo-Hookean, Mooney-Rivlin, Yeoh, and Anssari-Benam models are simultaneously calibrated using uniaxial tension and pure shear data. 
Plots compare the calibrated constitutive models with experimental data points.
"""

# =========================
#        Modules
# =========================

import tkinter as tk
import csv
import numpy as np
from tkinter import filedialog, messagebox, font
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.optimize import least_squares
try:
    import scienceplots
    plt.style.use(["science", "no-latex"])
except Exception:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 11,
        "axes.grid": True,
        "grid.linestyle": "--",
        "grid.alpha": 0.5,
    })


# =========================
#   Module-level helpers
# =========================

def _metrics(y_exp, y_model):
    res = y_exp - y_model
    rmse = np.sqrt(np.mean(res**2))
    nrmse = rmse / (y_exp.max() - y_exp.min()) * 100
    ss_res = np.sum(res**2)
    ss_tot = np.sum((y_exp - y_exp.mean())**2)
    r2 = 1 - ss_res / ss_tot
    max_err = np.max(np.abs(res))
    return rmse, nrmse, r2, max_err


def _build_rows(stress_u, stress_s, pred_u, pred_s):
    stress_comb = np.concatenate([stress_u, stress_s])
    pred_comb   = np.concatenate([pred_u,   pred_s])
    rows = []
    for label, y_e, y_m in [("Uniaxial",   stress_u,    pred_u),
                             ("Pure shear", stress_s,    pred_s),
                             ("Combined",   stress_comb, pred_comb)]:
        rows.append((label, *_metrics(y_e, y_m)))
    return rows


def _format_table(title, params_line, rows):
    header = f"{'Test mode':<12}{'RMSE':>12}{'NRMSE %':>12}{'R2':>10}{'Max error':>14}"
    sep = "-" * len(header)
    lines = [title, params_line, sep, header, sep]
    for label, rmse, nrmse, r2, max_err in rows:
        lines.append(f"{label:<12}{rmse:>12.4f}{nrmse:>12.2f}{r2:>10.4f}{max_err:>14.4f}")
    lines.append(sep)
    return "\n".join(lines)


# =========================
#      App definition
# =========================

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("My App")
        self.root.geometry("750x750")

        self.page1 = tk.Frame(root)
        self.page2 = tk.Frame(root)
        self.page3 = tk.Frame(root)
        self.page4 = tk.Frame(root)

        self.uniaxial_stretch = []
        self.uniaxial_stress  = []
        self.shear_stretch    = []
        self.shear_stress     = []

        self.selected_file  = None
        self.selected_model = None
        self.fit_params     = {}

        self.build_page1()
        self.build_page2()
        self.build_page3()

        self.page1.pack(fill="both", expand=True)

    # ------------------------------------------------------------------
    #  Pages 1–3
    # ------------------------------------------------------------------

    def build_page1(self):
        import webbrowser

        content_frame = tk.Frame(self.page1)
        content_frame.pack(fill="both", expand=True, padx=40, pady=35)

        text_widget = tk.Text(
            content_frame,
            wrap="word",
            width=60,
            height=14,
            borderwidth=0,
            highlightthickness=0,
            relief="flat"
        )

        bold_font = font.Font(text_widget, text_widget.cget("font"))
        bold_font.configure(weight="bold")

        text_widget.tag_configure("center",    justify="center")
        text_widget.tag_configure("bold_text", font=bold_font)
        text_widget.tag_configure("link",      foreground="blue", underline=True)
        text_widget.tag_bind("link", "<Button-1>", lambda e: webbrowser.open("https://mohsafaei.github.io"))
        text_widget.tag_bind("link", "<Enter>",    lambda e: text_widget.config(cursor="hand2"))
        text_widget.tag_bind("link", "<Leave>",    lambda e: text_widget.config(cursor=""))

        content_part1 = (
            "This app is designed for calibration of hyperelastic materials.\n"
            "You need to add your stress-strain data in the first place.\n"
            "Then, the app fits the model to the assigned data and yields the corresponding material parameters.\n"
            "Finally, the error metrics are yielded.\n\n"
            "-------------------------------------------\n\n"
            "The app is designed by\n"
        )
        my_name      = "Mohammad Ali Safaei"
        content_part2 = ".\nYou could check\n"
        website       = "mohsafaei.github.io"

        text_widget.insert("end", content_part1,  "center")
        text_widget.insert("end", my_name,         ("bold_text", "center"))
        text_widget.insert("end", content_part2,   "center")
        text_widget.insert("end", website,         ("bold_text", "link", "center"))

        text_widget.configure(state="disabled")
        text_widget.pack(pady=10)

        nav_frame = tk.Frame(self.page1)
        nav_frame.pack(side="bottom", fill="x", pady=15)
        tk.Button(nav_frame, text="Start", command=self.show_page2, width=20).pack()

    def browse_file(self):
        filepath = filedialog.askopenfilename(
            title="Select a file",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if not filepath:
            return False

        self.selected_file = filepath
        print(filepath)

        self.uniaxial_stretch = []
        self.uniaxial_stress  = []
        self.shear_stretch    = []
        self.shear_stress     = []

        try:
            with open(self.selected_file, "r", newline="", encoding="utf-8-sig") as csv_file:
                reader = csv.reader(csv_file)
                next(reader, None)

                for row_number, row in enumerate(reader, start=2):
                    if not row or all(not value.strip() for value in row):
                        continue

                    if len(row) < 4:
                        raise ValueError(f"Row {row_number} has fewer than four columns.")

                    self.uniaxial_stretch.append(float(row[0]))
                    self.uniaxial_stress.append(float(row[1]))
                    self.shear_stretch.append(float(row[2]))
                    self.shear_stress.append(float(row[3]))

            if not self.uniaxial_stretch:
                raise ValueError("The CSV file contains no data rows.")

        except (OSError, ValueError) as error:
            messagebox.showerror("Invalid CSV file", f"Could not read the file:\n{error}")
            return False

        self.show_page3()
        return True

    def build_page2(self):
        content_frame = tk.Frame(self.page2)
        content_frame.pack(fill="both", expand=True, padx=40, pady=30)

        label = tk.Label(
            content_frame,
            text=(
                "You need to browse your input data.\n"
                "the put file should be a .csv file containing the\n"
                "stress-strain results of uniaxial tension and pure shear tests.\n"
                "Please ensure that the columns in the selected CSV file follow this order:\n\n"
                "Column 1  →  Uniaxial stretch\n"
                "Column 2  →  Uniaxial stress\n"
                "Column 3  →  Shear stretch\n"
                "Column 4  →  Shear stress"
            ),
            wraplength=400,
            justify="center"
        )
        label.pack(pady=30)

        nav_frame = tk.Frame(self.page2)
        nav_frame.pack(side="bottom", fill="x", pady=15)
        tk.Button(nav_frame, text="← Back",  command=self.show_page1, width=12).pack(side="left",  padx=20)
        tk.Button(nav_frame, text="Browse",  command=self.browse_file, width=15).pack(side="right", padx=20)

    def build_page3(self):
        content_frame = tk.Frame(self.page3)
        content_frame.pack(fill="both", expand=True, padx=40, pady=20)

        label = tk.Label(
            content_frame,
            text=(
                "You need to choose your constitutive model:\n"
                "Neo-Hookean, Mooney-Rivlin, Yeoh, Gent and "
                "Anssari-Benam models are provided."
            ),
            wraplength=400,
            justify="center"
        )
        label.pack(pady=20)

        button_frame = tk.Frame(content_frame)
        button_frame.pack(pady=10)

        btn1 = tk.Button(button_frame, text="Neo-Hookean",   width=20, command=lambda: self.select_model("Neo-Hookean"))
        btn2 = tk.Button(button_frame, text="Mooney-Rivlin", width=20, command=lambda: self.select_model("Mooney-Rivlin"))
        btn3 = tk.Button(button_frame, text="Gent",          width=20, command=lambda: self.select_model("Gent"))
        btn4 = tk.Button(button_frame, text="Yeoh",          width=20, command=lambda: self.select_model("Yeoh"))
        btn5 = tk.Button(button_frame, text="Anssari-Benam", width=20, command=lambda: self.select_model("Anssari-Benam"))

        btn1.grid(row=0, column=0, padx=20, pady=8)
        btn2.grid(row=0, column=1, padx=20, pady=8)
        btn3.grid(row=1, column=0, padx=20, pady=8)
        btn4.grid(row=1, column=1, padx=20, pady=8)
        btn5.grid(row=2, column=0, padx=20, pady=8)

        nav_frame = tk.Frame(self.page3)
        nav_frame.pack(side="bottom", fill="x", pady=15)
        tk.Button(nav_frame, text="← Back", command=self.show_page2, width=12).pack(side="left", padx=20)

    # ------------------------------------------------------------------
    #  Model fitting
    # ------------------------------------------------------------------

    def select_model(self, model_name):
        self.selected_model = model_name
        print("Selected model:", model_name)

        lam_u    = np.array(self.uniaxial_stretch, dtype=float)
        stress_u = np.array(self.uniaxial_stress,  dtype=float)
        lam_s    = np.array(self.shear_stretch,    dtype=float)
        stress_s = np.array(self.shear_stress,     dtype=float)

        # ---------- Neo-Hookean ----------
        if model_name == "Neo-Hookean":
            def nh_uniaxial(lam, c10):
                return 2 * c10 * (lam - lam**(-2))

            def nh_pureshear(lam, c10):
                return 2 * c10 * (lam - lam**(-3))

            def residuals(x):
                c10 = x[0]
                r_u = stress_u[:8] - nh_uniaxial(lam_u[:8], c10)
                r_s = stress_s      - nh_pureshear(lam_s,    c10)
                return np.concatenate([r_u, r_s])

            result = least_squares(residuals, x0=[0.1], bounds=([1e-6], [5]))
            c10 = result.x[0]
            print(f"Neo-Hookean c10: {c10:.6g} MPa")
            self.fit_params = {"c10": c10}

        # ---------- Mooney-Rivlin ----------
        elif model_name == "Mooney-Rivlin":
            def mr_uniaxial(lam, c10, c01):
                return 2 * c10 * (lam - 1/lam**2) + 2 * c01 * (1 - 1/lam**3)

            def mr_pureshear(lam, c10, c01):
                return 2 * (lam - 1/lam**3) * (c10 + c01)

            def residuals(x):
                c10, c01 = x
                r_u = stress_u[:8] - mr_uniaxial(lam_u[:8], c10, c01)
                r_s = stress_s      - mr_pureshear(lam_s,    c10, c01)
                return np.concatenate([r_u, r_s])

            result = least_squares(residuals, x0=[0.1, 0.05], bounds=([1e-6, 1e-3], [5, 5]))
            c10, c01 = result.x
            print(f"Mooney-Rivlin c10, c01: {c10:.6g} MPa, {c01:.6g} MPa")
            self.fit_params = {"c10": c10, "c01": c01}

        # ---------- Yeoh ----------
        elif model_name == "Yeoh":
            def yeoh_uniaxial(lam, c1, c2, c3):
                i1  = lam**2 + 2/lam
                dW  = c1 + 2*c2*(i1-3) + 3*c3*(i1-3)**2
                return 2 * (lam - 1/lam**2) * dW

            def yeoh_pureshear(lam, c1, c2, c3):
                i1  = lam**2 + 1/lam**2 + 1
                dW  = c1 + 2*c2*(i1-3) + 3*c3*(i1-3)**2
                return 2 * (lam - 1/lam**3) * dW

            def residuals(x):
                c1, c2, c3 = x
                r_u = stress_u - yeoh_uniaxial(lam_u, c1, c2, c3)
                r_s = stress_s - yeoh_pureshear(lam_s, c1, c2, c3)
                return np.concatenate([r_u, r_s])

            result = least_squares(residuals, x0=[0.2, -0.015, 0.003],
                                   bounds=([-1, -1, -1], [1, 1, 1]))
            c1, c2, c3 = result.x
            print(f"Yeoh c1, c2, c3: {c1:.6g} MPa, {c2:.6g} MPa, {c3:.6g} MPa")
            self.fit_params = {"c1": c1, "c2": c2, "c3": c3}

        # ---------- Gent ----------
        elif model_name == "Gent":
            def gent_uniaxial(lam, mu, Jm):
                I1 = lam**2 + 2/lam
                return 2 * (lam - 1/lam**2) * mu / (2 * (1 - (I1-3)/Jm))

            def gent_pureshear(lam, mu, Jm):
                I1 = lam**2 + 1 + 1/lam**2
                return 2 * (lam - 1/lam**3) * mu / (2 * (1 - (I1-3)/Jm))

            def residuals(x):
                mu, Jm = x
                r_u = stress_u - gent_uniaxial(lam_u, mu, Jm)
                r_s = stress_s - gent_pureshear(lam_s, mu, Jm)
                return np.concatenate([r_u, r_s])

            result = least_squares(residuals, x0=[0.1, 10.0],
                                   bounds=([1e-6, 1e-3], [10.0, 500.0]))
            mu, Jm = result.x
            print(f"Gent μ: {mu:.6g} MPa, Jm: {Jm:.6g}")
            self.fit_params = {"mu": mu, "Jm": Jm}

        # ---------- Anssari-Benam ----------
        elif model_name == "Anssari-Benam":
            def dW_dI1(I1, mu1, N1, n1, beta1):
                A = 3*(n1-1)*mu1*N1 / (2*n1)
                return A * beta1 * ((I1-3)**(beta1-1) / (3*N1*(n1-1)) - 1/(I1-3*N1))

            def dW_dI2(I2, C20, eps1):
                return C20 * eps1 / 3 * (I2/3)**(eps1-1)

            def ab_uniaxial(lam, mu1, N1, n1, beta1, C20, eps1):
                I1 = lam**2 + 2/lam
                I2 = 2*lam  + 1/lam**2
                w1 = dW_dI1(I1, mu1, N1, n1, beta1)
                w2 = dW_dI2(I2, C20, eps1)
                return 2*(lam - 1/lam**2)*w1 + 2*(1 - 1/lam**3)*w2

            def ab_pureshear(lam, mu1, N1, n1, beta1, C20, eps1):
                I1 = lam**2 + 1 + 1/lam**2
                I2 = I1
                w1 = dW_dI1(I1, mu1, N1, n1, beta1)
                w2 = dW_dI2(I2, C20, eps1)
                return 2*(lam - 1/lam**3)*(w1 + w2)

            def residuals(x):
                mu1, N1, n1, beta1, C20, eps1 = x
                r_u = stress_u - ab_uniaxial(lam_u, mu1, N1, n1, beta1, C20, eps1)
                r_s = stress_s - ab_pureshear(lam_s, mu1, N1, n1, beta1, C20, eps1)
                return np.concatenate([r_u, r_s])

            result = least_squares(
                residuals,
                x0=[0.05, 0.9, 0.7, 2.5, 2.5, 0.25],
                bounds=(
                    [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6],
                    [10.0, 1.0,  10.0, 5.0,  10.0, 5.0 ]
                )
            )
            mu1, N1, n1, beta1, C20, eps1 = result.x
            print(f"Anssari-Benam μ, N, n, β, C₂₀, ε: "
                  f"{mu1:.4g} MPa, {N1:.7g}, {n1:.4g}, "
                  f"{beta1:.4g}, {C20:.4g} MPa, {eps1:.4g}")
            self.fit_params = {"mu1": mu1, "N1": N1, "n1": n1,
                               "beta1": beta1, "C20": C20, "eps1": eps1}

        else:
            messagebox.showinfo("Model", f"{model_name} is not implemented yet.")
            return

        self.save_metrics_table()
        self.build_page4()
        self.show_page4()

    # ------------------------------------------------------------------
    #  Metrics reporting  (real class method, uses self.fit_params)
    # ------------------------------------------------------------------

    def save_metrics_table(self):
        model    = self.selected_model
        lam_u    = np.array(self.uniaxial_stretch, dtype=float)
        stress_u = np.array(self.uniaxial_stress,  dtype=float)
        lam_s    = np.array(self.shear_stretch,    dtype=float)
        stress_s = np.array(self.shear_stress,     dtype=float)

        if model == "Neo-Hookean":
            c10    = self.fit_params["c10"]
            pred_u = 2 * c10 * (lam_u - lam_u**(-2))
            pred_s = 2 * c10 * (lam_s - lam_s**(-3))
            title       = "=== Neo-Hookean: joint calibration (uniaxial + pure shear) ==="
            params_line = f"c10 = {c10:.6g} MPa"

        elif model == "Mooney-Rivlin":
            c10, c01 = self.fit_params["c10"], self.fit_params["c01"]
            pred_u   = 2*c10*(lam_u - 1/lam_u**2) + 2*c01*(1 - 1/lam_u**3)
            pred_s   = 2*(lam_s - 1/lam_s**3) * (c10 + c01)
            title       = "=== Mooney-Rivlin: joint calibration (uniaxial + pure shear) ==="
            params_line = f"c10 = {c10:.6g} MPa,  c01 = {c01:.6g} MPa"

        elif model == "Yeoh":
            c1, c2, c3 = self.fit_params["c1"], self.fit_params["c2"], self.fit_params["c3"]
            i1_u   = lam_u**2 + 2/lam_u
            i1_s   = lam_s**2 + 1/lam_s**2 + 1
            dW_u   = c1 + 2*c2*(i1_u-3) + 3*c3*(i1_u-3)**2
            dW_s   = c1 + 2*c2*(i1_s-3) + 3*c3*(i1_s-3)**2
            pred_u = 2*(lam_u - 1/lam_u**2)*dW_u
            pred_s = 2*(lam_s - 1/lam_s**3)*dW_s
            title       = "=== Yeoh: joint calibration (uniaxial + pure shear) ==="
            params_line = f"c1 = {c1:.6g} MPa,  c2 = {c2:.6g} MPa,  c3 = {c3:.6g} MPa"

        elif model == "Gent":
            mu, Jm = self.fit_params["mu"], self.fit_params["Jm"]
            I1_u   = lam_u**2 + 2/lam_u
            I1_s   = lam_s**2 + 1 + 1/lam_s**2
            pred_u = 2*(lam_u - 1/lam_u**2) * mu / (2*(1-(I1_u-3)/Jm))
            pred_s = 2*(lam_s - 1/lam_s**3) * mu / (2*(1-(I1_s-3)/Jm))
            title       = "=== Gent: joint calibration (uniaxial + pure shear) ==="
            params_line = f"mu = {mu:.6g} MPa,  Jm = {Jm:.6g}"

        elif model == "Anssari-Benam":
            mu1   = self.fit_params["mu1"]
            N1    = self.fit_params["N1"]
            n1    = self.fit_params["n1"]
            beta1 = self.fit_params["beta1"]
            C20   = self.fit_params["C20"]
            eps1  = self.fit_params["eps1"]
            A     = 3*(n1-1)*mu1*N1 / (2*n1)

            I1_u  = lam_u**2 + 2/lam_u
            I2_u  = 2*lam_u  + 1/lam_u**2
            w1_u  = A*beta1*((I1_u-3)**(beta1-1)/(3*N1*(n1-1)) - 1/(I1_u-3*N1))
            w2_u  = C20*eps1/3 * (I2_u/3)**(eps1-1)
            pred_u = 2*(lam_u - 1/lam_u**2)*w1_u + 2*(1 - 1/lam_u**3)*w2_u

            I1_s  = lam_s**2 + 1 + 1/lam_s**2
            I2_s  = I1_s
            w1_s  = A*beta1*((I1_s-3)**(beta1-1)/(3*N1*(n1-1)) - 1/(I1_s-3*N1))
            w2_s  = C20*eps1/3 * (I2_s/3)**(eps1-1)
            pred_s = 2*(lam_s - 1/lam_s**3)*(w1_s + w2_s)

            title       = "=== Anssari-Benam: joint calibration (uniaxial + pure shear) ==="
            params_line = (f"mu1 = {mu1:.6g} MPa, N1 = {N1:.6g}, n1 = {n1:.6g}, "
                           f"beta1 = {beta1:.6g}, C20 = {C20:.6g} MPa, eps1 = {eps1:.6g}")

        else:
            return

        # metrics are computed against the full arrays (no [:8] slice here —
        # fitting used [:8] for NH/MR for numerical reasons, but metrics
        # should match the arrays used for plotting; adjust if you prefer [:8])
        rows  = _build_rows(stress_u, stress_s, pred_u, pred_s)
        table = _format_table(title, params_line, rows)
        print("\n" + table)

        results_path = f"{model}_fitting_result.txt"
        with open(results_path, "w") as f:
            f.write("Hyperelastic material calibration - fit results\n")
            f.write("=" * 60 + "\n\n")
            f.write(table + "\n")
        print(f"Saved: {results_path}")

    # ------------------------------------------------------------------
    #  Results page
    # ------------------------------------------------------------------

    def build_page4(self):
        for widget in self.page4.winfo_children():
            widget.destroy()

        fig = Figure(figsize=(8, 4), dpi=100)
        ax1, ax2 = fig.subplots(1, 2)

        lam_u    = np.array(self.uniaxial_stretch, dtype=float)
        stress_u = np.array(self.uniaxial_stress,  dtype=float)
        lam_s    = np.array(self.shear_stretch,    dtype=float)
        stress_s = np.array(self.shear_stress,     dtype=float)

        if self.selected_model == "Neo-Hookean":
            c10    = self.fit_params["c10"]
            pred_u = 2 * c10 * (lam_u - 1/lam_u**2)
            pred_s = 2 * c10 * (lam_s - 1/lam_s**3)
            label_text = "Neo-Hookean fit"

        elif self.selected_model == "Mooney-Rivlin":
            c10, c01 = self.fit_params["c10"], self.fit_params["c01"]
            pred_u   = 2*c10*(lam_u - 1/lam_u**2) + 2*c01*(1 - 1/lam_u**3)
            pred_s   = 2*(lam_s - 1/lam_s**3) * (c10 + c01)
            label_text = "Mooney-Rivlin fit"

        elif self.selected_model == "Yeoh":
            c1, c2, c3 = self.fit_params["c1"], self.fit_params["c2"], self.fit_params["c3"]
            i1_u   = lam_u**2 + 2/lam_u
            i1_s   = lam_s**2 + 1/lam_s**2 + 1
            dW_u   = c1 + 2*c2*(i1_u-3) + 3*c3*(i1_u-3)**2
            dW_s   = c1 + 2*c2*(i1_s-3) + 3*c3*(i1_s-3)**2
            pred_u = 2*(lam_u - 1/lam_u**2)*dW_u
            pred_s = 2*(lam_s - 1/lam_s**3)*dW_s
            label_text = "Yeoh fit"

        elif self.selected_model == "Gent":
            mu, Jm = self.fit_params["mu"], self.fit_params["Jm"]
            I1_u   = lam_u**2 + 2/lam_u
            I1_s   = lam_s**2 + 1 + 1/lam_s**2
            pred_u = 2*(lam_u - 1/lam_u**2) * mu / (2*(1-(I1_u-3)/Jm))
            pred_s = 2*(lam_s - 1/lam_s**3) * mu / (2*(1-(I1_s-3)/Jm))
            label_text = "Gent fit"

        elif self.selected_model == "Anssari-Benam":
            mu1   = self.fit_params["mu1"]
            N1    = self.fit_params["N1"]
            n1    = self.fit_params["n1"]
            beta1 = self.fit_params["beta1"]
            C20   = self.fit_params["C20"]
            eps1  = self.fit_params["eps1"]
            A     = 3*(n1-1)*mu1*N1 / (2*n1)

            I1_u  = lam_u**2 + 2/lam_u
            I2_u  = 2*lam_u  + 1/lam_u**2
            w1_u  = A*beta1*((I1_u-3)**(beta1-1)/(3*N1*(n1-1)) - 1/(I1_u-3*N1))
            w2_u  = C20*eps1/3 * (I2_u/3)**(eps1-1)
            pred_u = 2*(lam_u - 1/lam_u**2)*w1_u + 2*(1 - 1/lam_u**3)*w2_u

            I1_s  = lam_s**2 + 1 + 1/lam_s**2
            I2_s  = I1_s
            w1_s  = A*beta1*((I1_s-3)**(beta1-1)/(3*N1*(n1-1)) - 1/(I1_s-3*N1))
            w2_s  = C20*eps1/3 * (I2_s/3)**(eps1-1)
            pred_s = 2*(lam_s - 1/lam_s**3)*(w1_s + w2_s)
            label_text = "Anssari-Benam fit"

        else:
            messagebox.showerror("Error", "No valid model results are available.")
            return

        label_font = {"fontsize": 14, "fontweight": "bold"}
        title_font = {"fontsize": 15, "fontweight": "bold"}

        ax1.scatter(lam_u, stress_u, color="red", label="Data")
        ax1.plot(lam_u, pred_u, "b-", label=label_text)
        ax1.set_title("Uniaxial",   fontdict=title_font)
        ax1.set_xlabel("Stretch",   fontdict=label_font)
        ax1.set_ylabel("Stress",    fontdict=label_font)
        ax1.grid(True, alpha=0.2)
        ax1.legend()

        ax2.scatter(lam_s, stress_s, color="red", label="Data")
        ax2.plot(lam_s, pred_s, "b-", label=label_text)
        ax2.set_title("Pure Shear", fontdict=title_font)
        ax2.set_xlabel("Stretch",   fontdict=label_font)
        ax2.set_ylabel("Stress",    fontdict=label_font)
        ax2.grid(True, alpha=0.2)
        ax2.legend()

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.page4)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=(10, 5))

        parameter_lines = "\n".join(
            f"{k} = {v:.6g} MPa" for k, v in self.fit_params.items()
        )
        results_label = tk.Label(
            self.page4,
            text=f"Selected model: {self.selected_model}\nFitted material parameters:\n{parameter_lines}",
            justify="center",
            font=("Arial", 12, "bold")
        )
        results_label.pack(pady=(5, 5))

        nav_frame = tk.Frame(self.page4)
        nav_frame.pack(side="bottom", fill="x", pady=10)
        tk.Button(nav_frame, text="← Back",       command=self.show_page3, width=12).pack(side="left",  padx=20)
        tk.Button(nav_frame, text="Change Model",  command=self.show_page3, width=15).pack(side="right", padx=20)
        tk.Button(nav_frame, text="New File",      command=self.show_page2, width=12).pack(side="right", padx=5)

    # ------------------------------------------------------------------
    #  Navigation
    # ------------------------------------------------------------------

    def show_page1(self):
        self.page1.pack(fill="both", expand=True)
        self.page2.pack_forget()
        self.page3.pack_forget()
        self.page4.pack_forget()

    def show_page2(self):
        self.page1.pack_forget()
        self.page2.pack(fill="both", expand=True)
        self.page3.pack_forget()
        self.page4.pack_forget()

    def show_page3(self):
        self.page1.pack_forget()
        self.page2.pack_forget()
        self.page3.pack(fill="both", expand=True)
        self.page4.pack_forget()

    def show_page4(self):
        self.page1.pack_forget()
        self.page2.pack_forget()
        self.page3.pack_forget()
        self.page4.pack(fill="both", expand=True)


if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()

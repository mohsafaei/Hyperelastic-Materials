#!/usr/bin/env python3
"""grabit.py - a small, class-based image graph digitizer.
Requires Pillow; numpy and scipy are optional.
"""
import csv
import os
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk

try:
    from PIL import Image, ImageTk
except ImportError:  # give a useful error when launched
    Image = ImageTk = None
try:
    import numpy as np
except ImportError:
    np = None
try:
    from scipy.io import savemat
except ImportError:
    savemat = None


class GrabitApp:
    def __init__(self, root, filename=None):
        self.root = root
        root.title("grabit - Image Digitizer")
        root.geometry("1050x700")
        self.image = None
        self.tk_image = None
        self.image_path = None
        self.zoom = 1.0
        self.original_xlim = None
        self.original_ylim = None
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.mode = 'idle'
        self.pan_mode = False
        self.calib_pixel_pts = []       # X0, Xm, Y0, Ym, in image coordinates
        self.calib_values = []
        self.calib_step = 0
        self.x_log = tk.BooleanVar(value=False)
        self.y_log = tk.BooleanVar(value=False)
        self.im_dat = []
        self.true_dat = []
        self.datasets = {}
        self.pan_start = None
        self._build_ui()
        self.root.bind('<Key>', self.on_key_press)
        root.bind('<Configure>', lambda e: self.redraw())
        if filename:
            self.load_image(filename)

    def _build_ui(self):
        top = ttk.Frame(self.root, padding=4)
        top.pack(fill='x')
        ttk.Button(top, text='Load Image...', command=self.on_load_image).pack(side='left')
        ttk.Button(top, text='Restore Zoom', command=self.restore_zoom).pack(side='left', padx=3)
        ttk.Button(top, text='Calibrate', command=self.on_calibrate).pack(side='left', padx=3)
        ttk.Button(top, text='Cancel Calibration', command=self.cancel_calibration).pack(side='left')
        self.grab_button = ttk.Button(top, text='Grab Points', command=self.on_grab_points, state='disabled')
        self.grab_button.pack(side='left', padx=3)
        self.pan_button = ttk.Checkbutton(top, text='Pan', command=self.toggle_pan)
        self.pan_button.pack(side='left', padx=8)
        ttk.Checkbutton(top, text='X log', variable=self.x_log).pack(side='left')
        ttk.Checkbutton(top, text='Y log', variable=self.y_log).pack(side='left')
        ttk.Label(top, text='  a/b: zoom   Space: fit').pack(side='right')

        body = ttk.Frame(self.root)
        body.pack(fill='both', expand=True)
        self.canvas = tk.Canvas(body, background='#303030', highlightthickness=0, cursor='crosshair')
        self.canvas.pack(side='left', fill='both', expand=True)
        side = ttk.Frame(body, width=230, padding=6)
        side.pack(side='right', fill='y')
        ttk.Label(side, text='Datasets', font=('TkDefaultFont', 11, 'bold')).pack(anchor='w')
        self.listbox = tk.Listbox(side, height=18, exportselection=False)
        self.listbox.pack(fill='both', expand=True, pady=4)
        self.listbox.bind('<Double-Button-1>', self.on_dataset_double_click)
        buttons = ttk.Frame(side); buttons.pack(fill='x')
        ttk.Button(buttons, text='Save Current as Dataset', command=self.on_save_dataset).pack(fill='x')
        ttk.Button(buttons, text='Rename', command=self.on_rename_dataset).pack(fill='x', pady=2)
        ttk.Button(buttons, text='Delete', command=self.on_delete_dataset).pack(fill='x')
        ttk.Button(buttons, text='Export...', command=self.on_export_dataset).pack(fill='x', pady=2)
        self.status = ttk.Label(self.root, text='Idle', relief='sunken', anchor='w', padding=3)
        self.status.pack(side='bottom', fill='x')
        self.canvas.bind('<Button-1>', self.on_canvas_click)
        self.canvas.bind('<Button-3>', self.on_right_click)
        self.canvas.bind('<Button-2>', self.on_pan_start)
        self.canvas.bind('<B2-Motion>', self.on_pan_drag)
        self.canvas.bind('<ButtonRelease-2>', self.on_pan_end)
        self.canvas.bind('<MouseWheel>', self.on_mouse_wheel)
        self.canvas.bind('<Button-4>', lambda e: self.zoom_at(e, 1.15))
        self.canvas.bind('<Button-5>', lambda e: self.zoom_at(e, 1 / 1.15))

    def update_status(self, text):
        self.status.config(text=text)

    def on_load_image(self):
        if Image is None:
            messagebox.showerror('Missing dependency', 'Install Pillow: pip install Pillow')
            return
        p = filedialog.askopenfilename(filetypes=[('Images', '*.bmp *.jpg *.jpeg *.tif *.tiff *.gif *.png'), ('All files', '*.*')])
        if p:
            self.load_image(p)

    def load_image(self, path):
        try:
            self.image = Image.open(path).convert('RGB')
            self.image_path = path
            self.zoom = 1.0; self.pan_x = self.pan_y = 0
            self.cancel_calibration(); self.im_dat = []; self.true_dat = []
            self.fit_view()
            # Canvas equivalent of the original axis limits: the full image extent.
            self.original_xlim = (0, self.image.width)
            self.original_ylim = (0, self.image.height)
        except Exception as exc:
            messagebox.showerror('Load image', str(exc))

    def fit_view(self):
        if not self.image: return
        w, h = max(1, self.canvas.winfo_width()), max(1, self.canvas.winfo_height())
        self.zoom = min(w / self.image.width, h / self.image.height, 1.0)
        self.pan_x = (w - self.image.width * self.zoom) / 2
        self.pan_y = (h - self.image.height * self.zoom) / 2
        self.redraw()

    def redraw(self):
        self.canvas.delete('all')
        if not self.image: return
        w = max(1, int(self.image.width * self.zoom)); h = max(1, int(self.image.height * self.zoom))
        shown = self.image.resize((w, h), Image.Resampling.LANCZOS)
        self.tk_image = ImageTk.PhotoImage(shown)
        self.canvas.create_image(self.pan_x, self.pan_y, image=self.tk_image, anchor='nw', tags='image')
        self.redraw_markers()

    def image_to_canvas(self, p): return (self.pan_x + p[0] * self.zoom, self.pan_y + p[1] * self.zoom)
    def canvas_to_image(self, x, y): return ((x - self.pan_x) / self.zoom, (y - self.pan_y) / self.zoom)

    def redraw_markers(self):
        for i, p in enumerate(self.calib_pixel_pts):
            x, y = self.image_to_canvas(p); r = 5
            self.canvas.create_oval(x-r, y-r, x+r, y+r, outline='red', width=2, tags='marker')
            self.canvas.create_text(x+8, y-8, text=str(i+1), fill='red', anchor='sw', tags='marker')
        for i, p in enumerate(self.im_dat):
            x, y = self.image_to_canvas(p); r = 4
            self.canvas.create_oval(x-r, y-r, x+r, y+r, outline='red', fill='red', tags='point')
            if i: self.canvas.create_line(*self.image_to_canvas(self.im_dat[i-1]), x, y, fill='yellow', width=2, tags='point')

    def on_calibrate(self):
        if not self.image:
            messagebox.showinfo('Calibrate', 'Load an image first.'); return
        self.mode = 'calibration'; self.calib_pixel_pts = []; self.calib_values = []; self.calib_step = 0
        self.im_dat = []; self.true_dat = []; self.update_calibration_status(); self.redraw()

    def update_calibration_status(self):
        labels = ['ORIGIN of x-axis (X0)', 'MAXIMUM of x-axis (Xm)', 'ORIGIN of y-axis (Y0)', 'MAXIMUM of y-axis (Ym)']
        if self.mode == 'calibration': self.update_status('Calibration %d/4: Click on the %s' % (self.calib_step+1, labels[self.calib_step]))

    def cancel_calibration(self):
        self.mode = 'idle'; self.calib_step = 0; self.calib_pixel_pts = []; self.calib_values = []
        self.grab_button.config(state='normal' if self.is_calibrated() else 'disabled')
        self.update_status('Idle'); self.redraw()

    def on_canvas_click(self, event):
        if not self.image: return
        if self.pan_mode or self.mode == 'panning': return
        p = self.canvas_to_image(event.x, event.y)
        if self.mode == 'calibration':
            labels = ['Xo value', 'Xm value', 'Yo value', 'Ym value']
            value = simpledialog.askfloat(labels[self.calib_step], 'Enter real-world %s:' % labels[self.calib_step], parent=self.root)
            if value is None: return
            self.calib_pixel_pts.append(p); self.calib_values.append(value); self.calib_step += 1
            if self.calib_step == 4:
                if not self._validate_calibration():
                    self.calib_pixel_pts = []; self.calib_values = []; self.calib_step = 0; self.update_calibration_status(); return
                self.mode = 'idle'; self.grab_button.config(state='normal')
                self.update_status('Calibration complete. Click Grab Points.'); self.redraw()
            else: self.update_calibration_status(); self.redraw()
        elif self.mode == 'grabbing':
            self.im_dat.append(p)
            try: self.true_dat.append(self.pixel_to_data(*p))
            except ValueError as exc: messagebox.showerror('Transform', str(exc)); self.im_dat.pop(); return
            self.update_status('Grabbing points: %d points (Enter to finish)' % len(self.im_dat)); self.redraw_markers()

    def _validate_calibration(self):
        if self.calib_pixel_pts[0] == self.calib_pixel_pts[1] or self.calib_pixel_pts[2] == self.calib_pixel_pts[3]:
            messagebox.showerror('Calibration', 'Axis calibration points must be distinct.'); return False
        if self.x_log.get() and (self.calib_values[0] <= 0 or self.calib_values[1] <= 0):
            messagebox.showerror('Calibration', 'Log X values must be positive.'); return False
        if self.y_log.get() and (self.calib_values[2] <= 0 or self.calib_values[3] <= 0):
            messagebox.showerror('Calibration', 'Log Y values must be positive.'); return False
        return True

    def is_calibrated(self): return len(self.calib_pixel_pts) == 4 and len(self.calib_values) == 4

    def pixel_to_data(self, px, py):
        """Project a pixel onto each independently tilted axis, then interpolate.
        frac = dot(P-origin, axis_vector) / dot(axis_vector, axis_vector).
        This supports skew, rotation, mirrored axes, and log10 axis interpolation.
        """
        if not self.is_calibrated(): raise ValueError('Calibrate before transforming points.')
        x0, xm, y0, ym = self.calib_pixel_pts
        def frac(p, origin, end):
            vx, vy = end[0]-origin[0], end[1]-origin[1]
            return ((p[0]-origin[0])*vx + (p[1]-origin[1])*vy) / (vx*vx + vy*vy)
        fx = frac((px, py), x0, xm); fy = frac((px, py), y0, ym)
        def interp(f, a, b, logarithmic):
            if logarithmic:
                import math
                return 10 ** (math.log10(a) + f * (math.log10(b)-math.log10(a)))
            return a + f * (b-a)
        return (interp(fx, self.calib_values[0], self.calib_values[1], self.x_log.get()),
                interp(fy, self.calib_values[2], self.calib_values[3], self.y_log.get()))

    def on_grab_points(self):
        if not self.is_calibrated(): return
        self.mode = 'grabbing'; self.im_dat = []; self.true_dat = []
        self.update_status('Grabbing points: click points (Backspace/Delete removes last, Enter finishes)')
        self.canvas.focus_set(); self.redraw()

    def on_right_click(self, event):
        if self.mode != 'grabbing' or not self.im_dat: return
        p = self.canvas_to_image(event.x, event.y)
        i = min(range(len(self.im_dat)), key=lambda j: (self.im_dat[j][0]-p[0])**2 + (self.im_dat[j][1]-p[1])**2)
        self.im_dat.pop(i); self.true_dat.pop(i); self.redraw_markers()
        self.update_status('Grabbing points: %d points' % len(self.im_dat))

    def on_key_press(self, event):
        key = event.keysym.lower()
        if key in ('backspace', 'delete') and self.mode == 'grabbing' and self.im_dat:
            self.im_dat.pop(); self.true_dat.pop(); self.redraw_markers()
        elif key == 'return' and self.mode == 'grabbing':
            self.mode = 'idle'; self.update_status('Idle (points retained; save as dataset)'); self.redraw()
        elif key == 'a': self.zoom_at(self.canvas.winfo_width()/2, self.canvas.winfo_height()/2, 1.2)
        elif key == 'b': self.zoom_at(self.canvas.winfo_width()/2, self.canvas.winfo_height()/2, 1/1.2)
        elif key == 'space': self.fit_view()
        elif key == 'r': self.restore_zoom()

    def restore_zoom(self, event=None):
        """Reset the view back to the original loaded-image extent."""
        if self.original_xlim is None or self.original_ylim is None:
            self.update_status("No view to restore (load an image first)")
            return
        # This implementation uses a Tk Canvas rather than Matplotlib axes.
        self.fit_view()
        self.redraw_markers()
        self.canvas.draw_idle = getattr(self.canvas, 'draw_idle', lambda: None)
        self.canvas.draw_idle()
        self.update_status("View restored to original size")

    def zoom_at(self, event_or_x, factor, maybe_factor=None):
        if maybe_factor is None: x, y, f = self.canvas.winfo_width()/2, self.canvas.winfo_height()/2, factor
        else: x, y, f = event_or_x, factor, maybe_factor
        if hasattr(event_or_x, 'x'): x, y, f = event_or_x.x, event_or_x.y, factor
        if not self.image: return
        before = self.canvas_to_image(x, y); self.zoom = max(.05, min(20, self.zoom*f))
        self.pan_x = x - before[0]*self.zoom; self.pan_y = y - before[1]*self.zoom; self.redraw()

    def on_mouse_wheel(self, event): self.zoom_at(event, 1.15 if event.delta > 0 else 1/1.15)
    def toggle_pan(self):
        self.pan_mode = bool(self.pan_button.instate(['selected']))
        self.update_status('Panning (middle-drag or Pan mode)' if self.pan_mode else 'Idle')
    def on_pan_start(self, event): self.pan_start = (event.x, event.y, self.pan_x, self.pan_y); self.mode = 'panning'; self.update_status('Panning')
    def on_pan_drag(self, event):
        if self.pan_start:
            sx, sy, px, py = self.pan_start; self.pan_x = px+event.x-sx; self.pan_y = py+event.y-sy; self.redraw()
    def on_pan_end(self, event): self.pan_start = None; self.mode = 'idle'; self.update_status('Idle')

    def _selected_name(self):
        sel = self.listbox.curselection(); return self.listbox.get(sel[0]) if sel else None
    def on_save_dataset(self):
        if not self.true_dat: messagebox.showinfo('Dataset', 'There are no grabbed points.'); return
        name = simpledialog.askstring('Save dataset', 'Dataset name:', parent=self.root)
        if not name: return
        if name in self.datasets: messagebox.showerror('Dataset', 'That name already exists.'); return
        self.datasets[name] = [tuple(p) for p in self.true_dat]; self.listbox.insert('end', name)
    def on_rename_dataset(self):
        old = self._selected_name()
        if not old: return
        new = simpledialog.askstring('Rename dataset', 'New name:', initialvalue=old, parent=self.root)
        if not new or (new != old and new in self.datasets):
            if new and new in self.datasets: messagebox.showerror('Dataset', 'That name already exists.')
            return
        self.datasets[new] = self.datasets.pop(old); i = self.listbox.curselection()[0]; self.listbox.delete(i); self.listbox.insert(i, new); self.listbox.selection_set(i)
    def on_delete_dataset(self):
        name = self._selected_name()
        if name and messagebox.askyesno('Delete dataset', 'Delete %s?' % name):
            del self.datasets[name]; self.listbox.delete(self.listbox.curselection()[0])

    def on_export_dataset(self):
        name = self._selected_name()
        if not name: return
        types = [('Text', '*.txt'), ('CSV', '*.csv'), ('NumPy archive', '*.npz')]
        if savemat: types.append(('MATLAB MAT', '*.mat'))
        path = filedialog.asksaveasfilename(defaultextension='.txt', filetypes=types)
        if not path: return
        data = self.datasets[name]
        try:
            ext = os.path.splitext(path)[1].lower()
            if ext == '.csv':
                with open(path, 'w', newline='', encoding='utf-8') as f: csv.writer(f).writerows(data)
            elif ext == '.npz':
                if np is None: raise RuntimeError('NumPy is required for .npz export.')
                np.savez(path, **{'data': np.asarray(data, dtype=float)})
            elif ext == '.mat':
                if not savemat: raise RuntimeError('SciPy is not installed.')
                savemat(path, {'data': np.asarray(data, dtype=float) if np else data})
            else:
                with open(path, 'w', encoding='utf-8') as f:
                    for x, y in data: f.write('%s\t%s\n' % (x, y))
            self.update_status('Exported ' + os.path.basename(path))
        except Exception as exc: messagebox.showerror('Export', str(exc))

    def on_dataset_double_click(self, event=None):
        name = self._selected_name()
        if not name: return
        win = tk.Toplevel(self.root); win.title('Array Editor: '+name); win.geometry('420x350')
        tree = ttk.Treeview(win, columns=('x','y'), show='headings'); tree.heading('x', text='X'); tree.heading('y', text='Y'); tree.pack(fill='both', expand=True)
        for row in self.datasets[name]: tree.insert('', 'end', values=row)
        ttk.Button(win, text='Close', command=win.destroy).pack(pady=4)


def main():
    import sys
    if Image is None:
        print('Pillow is required: pip install Pillow'); return
    root = tk.Tk()
    app = GrabitApp(root, sys.argv[1] if len(sys.argv) > 1 else None)
    root.mainloop()


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Mesh + Contour viewer
- expects a folder `Output` in the same directory
- files in Output:
  - points : each line "r z" (two floats)
  - elements: each line "i1 i2 i3 i4" (four ints) - node indices: left-bottom, right-bottom, left-top, right-top
  - solution: each line "x y value" (three floats)

Run: python3 mesh_viewer.py
Dependencies: numpy, matplotlib
"""
import os
import sys
import numpy as np
import matplotlib

matplotlib.use('TkAgg')  # use Tk backend for embedding
from matplotlib.figure import Figure
from matplotlib.collections import PolyCollection
import matplotlib.tri as mtri
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk


# ---------------------- Utility functions ----------------------

def read_points(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            r = float(parts[0])
            z = float(parts[1])
            data.append((r, z))
    return np.array(data)  # shape (N,2)


def read_elements(path):
    elems = []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 4:
                continue
            elems.append([int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])])
    return np.array(elems, dtype=int)  # shape (M,4)


def read_solution(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            x = float(parts[0])
            y = float(parts[1])
            val = float(parts[2])
            data.append((x, y, val))
    return np.array(data)  # shape (K,3)


# ---------------------- Loading data ----------------------

def detect_and_load(output_dir='output'):
    pts_path = os.path.join(output_dir, 'points')
    elems_path = os.path.join(output_dir, 'elements')
    sol_path = os.path.join(output_dir, 'solution')

    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"Directory '{output_dir}' not found")
    for p in (pts_path, elems_path, sol_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"Required file '{p}' not found")

    pts = read_points(pts_path)
    elems = read_elements(elems_path)
    sol = read_solution(sol_path)

    # detect whether element indices are 1-based (common) or 0-based
    if elems.size == 0:
        raise ValueError('No elements read')
    if pts.size == 0:
        raise ValueError('No points read')

    if elems.min() == 1 and elems.max() <= pts.shape[0]:
        elems = elems - 1

    return pts, elems, sol


# ---------------------- Plotting helpers ----------------------

def build_polygons(pts, elems):
    # elems are [lb, rb, lt, rt] indexes
    polys = []
    for el in elems:
        # gather node coordinates in natural rectangular order
        lb = pts[el[0]]  # [x,y]
        rb = pts[el[1]]
        lt = pts[el[2]]
        rt = pts[el[3]]
        # polygon in either clockwise or ccw order
        poly = np.array([lb, rb, rt, lt])
        polys.append(poly)
    return np.array(polys)  # shape (M,4,2)


def calculate_element_centers(pts, elems):
    """Calculate center points for all elements"""
    centers = []
    for el in elems:
        # Get coordinates of all four nodes
        nodes = pts[el]
        # Calculate center as mean of all nodes
        center = np.mean(nodes, axis=0)
        centers.append(center)
    return np.array(centers)


def plot_mesh(ax, pts, elems, show_nodes=True, show_centers=True):
    polys = build_polygons(pts, elems)

    # Use PolyCollection which accepts list/array of (N,2) coordinate arrays
    pc = PolyCollection(polys, edgecolors='k', facecolors='none', linewidths=0.6)
    ax.add_collection(pc)

    # Plot nodes
    if show_nodes:
        ax.scatter(pts[:, 0], pts[:, 1], c='black', s=20, alpha=0.7)

    # Plot element centers
    if show_centers:
        centers = calculate_element_centers(pts, elems)
        ax.scatter(centers[:, 0], centers[:, 1], c='black', s=20, alpha=0.7)

    # set limits
    ax.autoscale_view()
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', labelsize=14)

    # Add legend if any points are shown
    if show_nodes or show_centers:
        ax.legend(loc='upper right', fontsize=10)


def plot_contour(ax, sol, app=None, show_nodes=True, show_centers=True, pts=None, elems=None):
    """Draw contour on ax from sol (Kx3 array). If `app` provided, use app.cbar to manage the colorbar
    so we don't create multiple colorbars and leave invisible axes that shift the layout.
    """
    # sol: (K,3) array -> x, y, value
    x = sol[:, 0]
    y = sol[:, 1]
    v = sol[:, 2]

    # Use triangulation-based contour (works for scattered data)
    tri = mtri.Triangulation(x, y)
    cf = ax.tricontourf(tri, v, cmap='jet', alpha=0.8)

    # Plot nodes if requested and data provided
    if show_nodes and pts is not None:
        ax.scatter(pts[:, 0], pts[:, 1], c='red', s=15, alpha=0.8,
                   label=f'Узлы ({len(pts)})', zorder=3)

    # Plot element centers if requested and data provided
    if show_centers and pts is not None and elems is not None:
        centers = calculate_element_centers(pts, elems)
        ax.scatter(centers[:, 0], centers[:, 1], c='white', s=12, alpha=0.8,
                   marker='s', edgecolors='black', linewidth=0.5,
                   label=f'Центры ({len(centers)})', zorder=3)

    # If an app object provided, remove previous colorbar safely
    if app is not None:
        try:
            if getattr(app, 'cbar', None) is not None:
                app.cbar.remove()
                app.cbar = None
        except Exception:
            # be defensive: if removing fails for any reason, ignore and continue
            app.cbar = None

    # create new colorbar and store it on app (if available) so it can be removed next time
    cb = ax.figure.colorbar(cf, ax=ax, orientation='vertical')
    if app is not None:
        app.cbar = cb

    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', labelsize=14)
    cb.ax.tick_params(labelsize=14)

    # Add legend if any points are shown
    if (show_nodes and pts is not None) or (show_centers and pts is not None and elems is not None):
        ax.legend(loc='upper right', fontsize=10)


class MeshApp:
    def __init__(self, master, pts, elems, sol):
        self.master = master
        self.pts = pts
        self.elems = elems
        self.sol = sol
        # reference for current colorbar (so we can remove it cleanly)
        self.cbar = None

        # Settings for displaying points
        self.show_nodes = tk.BooleanVar(value=True)
        self.show_centers = tk.BooleanVar(value=True)

        master.title('Mesh & Contour Viewer')

        # Main frame
        self.mainframe = ttk.Frame(master, padding=(2, 2, 2, 2))
        self.mainframe.pack(fill=tk.BOTH, expand=True)

        # Left: figure (large)
        self.fig = Figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.mainframe)
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Right: controls
        ctrl_frame = ttk.Frame(self.mainframe, width=200)
        ctrl_frame.pack(side=tk.RIGHT, fill=tk.Y)

        # View mode selection
        view_label = ttk.Label(ctrl_frame, text="Режим отображения:", font=('Arial', 10, 'bold'))
        view_label.pack(anchor=tk.NW, pady=(20, 5), padx=10)

        self.view_var = tk.StringVar(value='mesh')
        rb1 = ttk.Radiobutton(ctrl_frame, text='Сетка', variable=self.view_var, value='mesh',
                              command=self.on_mode_change)
        rb2 = ttk.Radiobutton(ctrl_frame, text='Контур', variable=self.view_var, value='contour',
                              command=self.on_mode_change)
        rb1.pack(anchor=tk.NW, pady=5, padx=10)
        rb2.pack(anchor=tk.NW, pady=5, padx=10)

        # Separator
        separator = ttk.Separator(ctrl_frame, orient='horizontal')
        separator.pack(fill='x', pady=10, padx=5)

        # Points display options
        points_label = ttk.Label(ctrl_frame, text="Отображение точек:", font=('Arial', 10, 'bold'))
        points_label.pack(anchor=tk.NW, pady=(10, 5), padx=10)

        cb1 = ttk.Checkbutton(ctrl_frame, text='Показывать узлы',
                              variable=self.show_nodes, command=self.on_display_change)
        cb2 = ttk.Checkbutton(ctrl_frame, text='Показывать центры',
                              variable=self.show_centers, command=self.on_display_change)
        cb1.pack(anchor=tk.NW, pady=2, padx=15)
        cb2.pack(anchor=tk.NW, pady=2, padx=15)

        # Info label
        self.info_label = ttk.Label(ctrl_frame, text="", font=('Arial', 9))
        self.info_label.pack(anchor=tk.NW, pady=(15, 5), padx=10)

        # initial draw
        self.update_info()
        self.draw()

    def clear_ax(self):
        self.ax.cla()

    def update_info(self):
        """Update information label"""
        nodes_count = len(self.pts)
        elements_count = len(self.elems)
        info_text = f"Узлы: {nodes_count}\nЭлементы: {elements_count}"
        self.info_label.config(text=info_text)

    def draw(self):
        mode = self.view_var.get()
        # remove any existing colorbar before clearing axes to avoid orphaned axes/layout shifts
        try:
            if getattr(self, 'cbar', None) is not None:
                self.cbar.remove()
                self.cbar = None
        except Exception:
            self.cbar = None

        self.clear_ax()

        if mode == 'mesh':
            plot_mesh(self.ax, self.pts, self.elems,
                      show_nodes=self.show_nodes.get(),
                      show_centers=self.show_centers.get())
        else:
            plot_contour(self.ax, self.sol, app=self,
                         show_nodes=self.show_nodes.get(),
                         show_centers=self.show_centers.get(),
                         pts=self.pts, elems=self.elems)

        self.ax.set_xlabel('R', fontsize=16)
        self.ax.set_ylabel('Z', fontsize=16)
        self.canvas.draw_idle()

    def on_mode_change(self):
        self.draw()

    def on_display_change(self):
        self.draw()


# ---------------------- Main ----------------------

def main():
    try:
        pts, elems, sol = detect_and_load('output')
    except Exception as e:
        print('Ошибка при загрузке данных:', e)
        sys.exit(1)

    root = tk.Tk()
    # make the window reasonably large
    root.geometry('1000x700')
    app = MeshApp(root, pts, elems, sol)
    root.mainloop()


if __name__ == '__main__':
    main()
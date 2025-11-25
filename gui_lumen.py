#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GUI that loads a FITS image, displays its dimensions,
lets the user choose how to divide the image using divisors
of height and width, and then processes selected blocks.

NEW FEATURES:
- Shows image dimensions
- Computes valid divisors for height and width
- User selects number of vertical and horizontal divisions
- Image is divided dynamically based on user choices
"""

import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.nddata import Cutout2D

from lumen_max import Lumen_max
from lumen_segment import lumen_segment_nobj
from lumen_segmentnobj1 import lumen_segment_nobj1


# ============================
# Utility: get all divisors
# ============================
def divisors(n):
    """Returns all positive divisors of n."""
    return [d for d in range(1, n + 1) if n % d == 0]


# ============================
# Processing function
# ============================
def procesar_bloque(img_data, block_index, blocks_y, blocks_x,
                    analysis_type, n_objects, fixed_radius=None):

    H, W = img_data.shape

    # Compute block size automatically
    block_h = H // blocks_y
    block_w = W // blocks_x

    # Compute row/col of block
    r = block_index // blocks_x
    c = block_index % blocks_x

    # Extract block
    bloque = img_data[r*block_h:(r+1)*block_h, c*block_w:(c+1)*block_w]

    # Central 400×400 crop for analysis
    cutout = Cutout2D(bloque, (block_w//2, block_h//2), (400, 400))
    cout = cutout.data

    imagen_final, maxima = Lumen_max(cout, n_objects)

    # Segmentation
    if analysis_type == "Estático":
        binary_mask, individual_masks = lumen_segment_nobj1(cout, n_objects, maxima, fixed_radius)
    else:
        binary_mask, individual_masks = lumen_segment_nobj(cout, n_objects, maxima, Radius='adaptive')

    return cout, maxima, binary_mask, individual_masks


# ============================
# GUI
# ============================
def cargar_fits():
    filepath = filedialog.askopenfilename(filetypes=[("FITS files", "*.fits")])
    if not filepath:
        return

    try:
        global img_data, H, W

        # Load image
        img_data = fits.getdata(filepath)
        H, W = img_data.shape

        # Update GUI with size
        status_var.set(f"Loaded image: {H} × {W}")

        # Compute divisors for user
        div_y = divisors(H)
        div_x = divisors(W)

        menu_y['values'] = div_y
        menu_x['values'] = div_x

        menu_y.current(0)
        menu_x.current(0)

        btn_setDivisions.config(state="normal")

    except Exception as e:
        messagebox.showerror("Error", str(e))


def establecer_divisiones():
    """Updates the block number dropdown once the user selects divisions."""
    global blocks_y, blocks_x

    blocks_y = int(menu_y.get())
    blocks_x = int(menu_x.get())
    blocks_total = blocks_y * blocks_x

    block_menu['values'] = list(range(blocks_total))
    block_menu.current(0)

    status_var.set(
        f"Image {H}×{W} divided into {blocks_y} × {blocks_x} blocks.\n"
        f"Block size = {H//blocks_y} × {W//blocks_x}"
    )

    btn_procesar.config(state="normal")


def ejecutar_proceso():
    try:
        block_index = int(block_menu.get())
        analysis_type = tipo_var.get()
        n_objects = int(entry_n.get())

        if analysis_type == "Estático":
            radius = int(entry_r.get())
        else:
            radius = None

        cout, maxima, binary_mask, individual_masks = procesar_bloque(
            img_data, block_index, blocks_y, blocks_x,
            analysis_type, n_objects, radius
        )

        # Plot results
        fig, axs = plt.subplots(1, 2, figsize=(8, 4))

        axs[0].imshow(cout, cmap='gray')
        axs[0].plot(maxima[:, 1], maxima[:, 0], 'r+', markersize=10)
        axs[0].set_title("Maxima")
        axs[0].axis('off')

        axs[1].imshow(binary_mask, cmap='gray')
        axs[1].set_title("Mask")
        axs[1].axis('off')

        canvas = FigureCanvasTkAgg(fig, master=frame_plot)
        canvas.draw()
        canvas.get_tk_widget().pack()

    except Exception as e:
        messagebox.showerror("Processing Error", str(e))


# ============================
# Build GUI
# ============================
root = tk.Tk()
root.title("FITS Block Selector - Divisor-based GUI")

frame_controls = tk.Frame(root)
frame_controls.pack(padx=10, pady=10, fill='x')

# Load FITS
btn_cargar = tk.Button(frame_controls, text="Load FITS Image", command=cargar_fits)
btn_cargar.grid(row=0, column=0, padx=5, pady=5)

status_var = tk.StringVar()
tk.Label(frame_controls, textvariable=status_var).grid(row=0, column=1, columnspan=4)

# Divisions selection
tk.Label(frame_controls, text="Vertical divisions:").grid(row=1, column=0)
menu_y = ttk.Combobox(frame_controls, width=5, state="readonly")
menu_y.grid(row=1, column=1)

tk.Label(frame_controls, text="Horizontal divisions:").grid(row=1, column=2)
menu_x = ttk.Combobox(frame_controls, width=5, state="readonly")
menu_x.grid(row=1, column=3)

btn_setDivisions = tk.Button(frame_controls, text="Set Divisions", command=establecer_divisiones, state='disabled')
btn_setDivisions.grid(row=2, column=0, columnspan=4, pady=5)


# Block selection
tk.Label(frame_controls, text="Block:").grid(row=3, column=0)
block_menu = ttk.Combobox(frame_controls, width=5, state="readonly")
block_menu.grid(row=3, column=1)

# Analysis parameters
tk.Label(frame_controls, text="Analysis:").grid(row=3, column=2)
tipo_var = tk.StringVar(value="Adaptativo")
ttk.Combobox(frame_controls, textvariable=tipo_var,
             values=["Adaptativo", "Estático"], state="readonly").grid(row=3, column=3)

tk.Label(frame_controls, text="Objects:").grid(row=4, column=0)
entry_n = tk.Entry(frame_controls, width=5)
entry_n.insert(0, "10")
entry_n.grid(row=4, column=1)

tk.Label(frame_controls, text="Radius (static):").grid(row=4, column=2)
entry_r = tk.Entry(frame_controls, width=5)
entry_r.insert(0, "10")
entry_r.grid(row=4, column=3)

btn_procesar = tk.Button(frame_controls, text="Process Block",
                         command=ejecutar_proceso, state='disabled')
btn_procesar.grid(row=5, column=0, columnspan=4, pady=10)

# Plot frame
frame_plot = tk.Frame(root)
frame_plot.pack(fill='both', expand=True)

root.mainloop()

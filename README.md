FITS Bright Object Detection GUI

This repository contains a Python-based graphical user interface designed for the detection and segmentation of bright objects within large astronomical FITS images. The tool enables users to load FITS files, visualize their structure, divide them into user-defined blocks, and apply adaptive or static segmentation methods to isolate luminous regions of interest.

The application is intended for research workflows where high-resolution images—such as 12,288 × 12,288-pixel mosaics—must be processed efficiently and interactively. By allowing the user to segment the image into manageable blocks, the tool significantly reduces processing time and memory usage.

Features

FITS Image Loading
Supports standard FITS images through astropy.io.fits.

Automatic Dimension Detection
Reads image dimensions on load and reports size, shape, and total pixel count.

User-Defined Block Division
The user may choose any block size that divides evenly into the image’s width and height.

Central Cutout Extraction
Each block is reduced to a 400 × 400 central region for focused analysis.

Maximum Detection Algorithm (Lumen_max)
Identifies local maxima within the cutout, representing bright astronomical sources.

Two Segmentation Modes

Adaptive Radius: Automatically adjusts object radius.

Static Radius: Uses a user-defined circular mask radius.

Interactive Visualization
The GUI displays maxima positions and their corresponding binary segmentation masks.

Backend Algorithms Included
Custom segmentation routines (lumen_segment_nobj, lumen_segment_nobj1) included in the repository.

Requirements

The project uses:

Python 3.8+

Tkinter (bundled with Python)

NumPy

Matplotlib

Astropy

Astropy Cutout2D utilities

tkinter.ttk widgets


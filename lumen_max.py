import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from skimage import exposure
import cv2

def Lumen_max(subs_img, n=10):
    """
    Detects and labels the brightest regional maxima in an image.

    Parameters
    ----------
    subs_img : 2D numpy array
        Input image (typically a cropped FITS section).
    n : int, optional
        Number of brightest peaks to extract. Default is 10.

    Returns
    -------
    thresholded : numpy array (H x W x 3, uint8)
        RGB visualization of the image with the top `n` brightest objects labeled.
    maxima : numpy array (n x 2)
        Array containing the (row, col) coordinates of the selected peaks.

    Notes
    -----
    - Local maxima are detected using a 3×3 footprint.
    - Peaks are sorted by intensity and filtered to enforce a minimum distance.
    - The top `n` brightest maxima are overlayed using OpenCV.
    - A 3-panel visualization is shown using Matplotlib.
    """

    # --- Local maxima detection ---
    footprint = np.ones((3, 3))              # Neighborhood for peak detection
    coordinates = peak_local_max(subs_img,   # Returns array of (row, col)
                                 footprint=footprint)

    # Boolean mask of detected regional maxima
    regional_max = np.zeros_like(subs_img, dtype=bool)
    regional_max[tuple(coordinates.T)] = True

    # Extract their coordinates and intensities
    rows, cols = np.where(regional_max)
    intensities = subs_img[rows, cols]

    # Sort maxima by brightness (descending)
    sort_idx = np.argsort(intensities)[::-1]
    sorted_intensities = intensities[sort_idx]
    sorted_positions = np.column_stack((rows[sort_idx], cols[sort_idx]))

    # --- Peak selection with minimum distance constraint ---
    min_dist = 10         # Minimum spacing between selected peaks
    maxima = []

    while len(sorted_positions) > 0 and len(maxima) < n:
        current_peak = sorted_positions[0]
        maxima.append(current_peak)

        # Remove peaks too close to the chosen one
        distances = np.linalg.norm(sorted_positions - current_peak, axis=1)
        keep_idx = distances > min_dist
        sorted_positions = sorted_positions[keep_idx]
        sorted_intensities = sorted_intensities[keep_idx]

    # If fewer than n found, pad with remaining peaks
    if len(maxima) < n:
        needed = n - len(maxima)
        maxima.extend(sorted_positions[:needed])

    maxima = np.array(maxima)

    # --- Prepare RGB visualization image ---
    thresholded = np.stack([subs_img]*3, axis=-1)  # Convert to 3-channel grayscale

    # Contrast stretching (1st to 99th percentile)
    thresholded = exposure.rescale_intensity(
        thresholded,
        in_range=(np.percentile(thresholded, 1), np.percentile(thresholded, 99))
    )
    thresholded = (thresholded * 255).astype(np.uint8)

    # --- Drawing parameters ---
    marker_size = max(3, int(min(subs_img.shape) / 60))
    text_offset = marker_size * 2
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = 0.2
    font_thickness = 1

    # Overlay markers and numeric labels
    for i, (row, col) in enumerate(maxima, start=1):

        center = (int(col), int(row))  # OpenCV uses (x, y)

        # Draw filled red circle with alpha blending
        overlay = thresholded.copy()
        cv2.circle(overlay, center, marker_size, (0, 0, 255), -1)
        alpha = 0.7
        cv2.addWeighted(overlay, alpha, thresholded, 1 - alpha, 0, thresholded)

        # Draw label background box and number
        text = str(i)
        ((text_w, text_h), _) = cv2.getTextSize(text, font, font_scale, font_thickness)
        text_origin = (center[0] + text_offset, center[1] - text_offset)

        cv2.rectangle(
            thresholded,
            (text_origin[0], text_origin[1] - text_h),
            (text_origin[0] + text_w, text_origin[1] + 5),
            (0, 0, 255),
            thickness=-1
        )
        cv2.putText(thresholded, text, (text_origin[0], text_origin[1]),
                    font, font_scale, (255, 255, 255), font_thickness)

    # --- Visualization (3-panel figure) ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    axes[0].imshow(subs_img, cmap='gray')
    axes[0].set_title('Processed Image')
    axes[0].axis('off')

    axes[1].imshow(subs_img, cmap='gray')
    axes[1].plot(coordinates[:, 1], coordinates[:, 0], 'r.', markersize=5)
    axes[1].set_title('All Local Maxima Detected')
    axes[1].axis('off')

    axes[2].imshow(cv2.cvtColor(thresholded, cv2.COLOR_BGR2RGB))
    axes[2].set_title(f'Top {n} Brightest Objects')
    axes[2].axis('off')

    plt.tight_layout()
    plt.show()

    return thresholded, maxima

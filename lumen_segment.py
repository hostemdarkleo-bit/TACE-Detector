import numpy as np
from scipy.ndimage import gaussian_filter1d, map_coordinates
from skimage.morphology import remove_small_objects
from scipy.ndimage import binary_fill_holes

def lumen_segment_nobj(img, n, maxima, Radius='adaptive', MinRadius=5, MaxRadius=50):
    """
    Generates circular binary masks for the n brightest detected objects.

    Parameters:
        img (2D np.array):
            Input image from which the objects were detected.

        n (int):
            Number of objects to segment.

        maxima (np.array):
            Array of coordinates [row, col] of the detected maxima.

        Radius (str or int):
            If 'adaptive', the radius is estimated per object.
            If numeric, a fixed radius is used for all objects.

        MinRadius (int):
            Minimum possible radius when using adaptive segmentation.

        MaxRadius (int):
            Maximum possible radius when using adaptive segmentation.

    Returns:
        binary_mask (2D bool np.array):
            Global mask containing all segmented objects.

        individual_masks (list of 2D bool np.array):
            One binary mask per segmented object.
    """

    # Initialize the global mask and list of individual masks
    binary_mask = np.zeros(img.shape, dtype=bool)
    individual_masks = []

    # Process the top n maxima (or fewer if maxima < n)
    for i in range(min(n, maxima.shape[0])):
        center = maxima[i]  # (row, col) of current object

        # Determine radius: fixed or adaptively estimated
        if Radius == 'adaptive':
            radius = estimate_optimal_radius(img, center, MinRadius, MaxRadius)
        else:
            radius = Radius

        # Generate a circular mask for the object
        yy, xx = np.indices(img.shape)
        mask = (xx - center[1])**2 + (yy - center[0])**2 <= radius**2

        individual_masks.append(mask)
        binary_mask = binary_mask | mask  # Add to combined mask

    # Fill holes inside segmented objects
    binary_mask = binary_fill_holes(binary_mask)

    # Remove undesired small regions or noise
    binary_mask = remove_small_objects(binary_mask, MinRadius**2)

    return binary_mask, individual_masks


def estimate_optimal_radius(img, center, minR, maxR):
    """
    Estimates an optimal radius for a bright region using radial intensity profiles.

    Parameters:
        img (2D np.array):
            Source image.

        center (tuple):
            (row, col) location of the object.

        minR (int):
            Minimum allowed radius.

        maxR (int):
            Maximum allowed radius.

    Returns:
        radius (int):
            Estimated optimal radius for the object.
    """

    # Sample radial directions (16 evenly spaced angles)
    theta = np.linspace(0, 2*np.pi, 16, endpoint=False)

    # Matrix to store intensity profiles for each angle
    intensity_profiles = np.zeros((len(theta), maxR))

    # Collect intensity values along each radial direction
    for i, angle in enumerate(theta):

        # Generate radial projection coordinates
        x_proj = center[1] + np.cos(angle) * np.arange(1, maxR + 1)
        y_proj = center[0] + np.sin(angle) * np.arange(1, maxR + 1)

        # Filter points that remain inside image boundaries
        valid = (
            (x_proj >= 0) & (x_proj < img.shape[1]) &
            (y_proj >= 0) & (y_proj < img.shape[0])
        )
        x_proj = x_proj[valid]
        y_proj = y_proj[valid]

        # Sample the image intensities along the projection
        intensity_profiles[i, :len(x_proj)] = map_coordinates(
            img, [y_proj, x_proj], order=1, mode='nearest'
        )

    # Average the intensity over all radial directions
    mean_profile = np.mean(intensity_profiles, axis=0)

    # Normalize the radial profile
    normalized_profile = mean_profile / np.max(mean_profile)

    # Smooth to reduce noise
    smoothed_profile = gaussian_filter1d(normalized_profile, sigma=1)

    # Threshold to detect intensity drop-off (edge of the object)
    threshold = np.percentile(smoothed_profile, 70)

    # Find first radius where intensity falls below threshold
    drop_points = np.where(smoothed_profile < threshold)[0]

    if drop_points.size == 0:
        # If no drop detected, use maximum radius
        radius = maxR
    else:
        # Clamp radius within allowed limits
        radius = max(minR, min(maxR, drop_points[0]))

    return radius

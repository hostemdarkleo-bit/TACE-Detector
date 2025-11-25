def lumen_segment_nobj1(img, n, maxima, Radius):
    """
    Generates circular binary masks for the n brightest objects,
    using a fixed radius for all objects.

    Parameters:
        img (2D np.array):
            Input image from which the objects were detected.

        n (int):
            Number of objects to segment.

        maxima (np.array):
            Array of coordinates [row, col] for the detected maxima.

        Radius (int):
            Fixed radius of the circular masks.

    Returns:
        binary_mask (2D bool np.array):
            Combined binary mask including all segmented objects.

        individual_masks (list of 2D bool np.array):
            One binary mask per segmented object.
    """

    import numpy as np

    # Create an empty global mask with same shape as the image
    binary_mask = np.zeros_like(img, dtype=bool)

    # List to store individual object masks
    individual_masks = []
    
    # Process each detected maximum up to n objects
    for i in range(min(n, maxima.shape[0])):
        center = maxima[i]  # (row, col) coordinates of current object
        
        # Create coordinate grids for distance calculation
        rr, cc = np.ogrid[:img.shape[0], :img.shape[1]]

        # Generate circular mask: points within the radius are True
        mask = (rr - center[0])**2 + (cc - center[1])**2 <= Radius**2
        
        # Store individual mask
        individual_masks.append(mask)

        # Update global mask by combining with OR operation
        binary_mask = binary_mask | mask
    
    return binary_mask, individual_masks

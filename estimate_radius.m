function radius = estimate_radius(img, center)
    % Radial intensity profile analysis
    [rows, cols] = size(img);
    max_radius = min([center(1)-1, rows-center(1), center(2)-1, cols-center(2)]);
    max_radius = min(max_radius, 50); % Limit maximum radius
    
    intensity_profile = zeros(1, max_radius);
    for r = 1:max_radius
        mask = create_mask(size(img), center, r);
        intensity_profile(r) = mean(img(mask));
    end
    
    % Find radius where intensity drops to 50% of maximum
    norm_profile = intensity_profile / max(intensity_profile);
    drop_point = find(norm_profile < 0.5, 1);
    
    if isempty(drop_point)
        radius = max_radius;
    else
        radius = max(5, min(drop_point, max_radius)); % Keep between 5 and max_radius
    end
end

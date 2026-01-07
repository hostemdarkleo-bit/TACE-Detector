function mask = create_mask(img_size, center, radius)
    [xx,yy] = meshgrid(1:img_size(2), 1:img_size(1));
    mask = sqrt((xx-center(2)).^2 + (yy-center(1)).^2) <= radius;
end

function [grain_size] = fn_calcGrainSizeImage(I2)
I2_vector = I2(:);
[unique_elements, ~, occurrence_count] = unique(I2_vector);
counts = accumarray(occurrence_count, 1);
area = counts;
grain_size = sqrt(area/pi)*2;
end
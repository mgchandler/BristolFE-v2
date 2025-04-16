function [mx, mn] = fn_estimate_max_min_vels_2(matls) %EDITED
v = [0,0];
j = 1;

    current_matl = matls{1}; %to solve indexing issue
    s = [];
    if isfield(current_matl, 'stiffness_matrix')
        s = current_matl.stiffness_matrix(:);
    end
    if isfield(current_matl, 'D') %choosing first matls
        s = current_matl.D(:);
    end
    if isfield(current_matl, 'density') %deal with legacy naming
        rho = current_matl.density;
    else
        rho = current_matl.rho;
        %rho = matls(i).rho;
    end

    if ~isempty(s)
        s = s(abs(s) > 0);
        v(j,1) = sqrt(min(abs(s)) / rho);
        v(j,2) = sqrt(max(abs(s)) / rho);
        j = j + 1;
    end



mn = min(v(:,1));
mx = max(v(:,2));
end
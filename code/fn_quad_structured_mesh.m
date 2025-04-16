function mod = fn_quad_structured_mesh(el_size,sample_width,sample_depth)
%SUMMARY
%   Utility function for generating a isometric structured mesh of triangular
%   elements, that fills the region specified by bdry_nds.
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 2 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 3 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%--------------------------------------------------------------------------

columns=ceil(sample_width/el_size); % ceil function rounds each element to nearest integer greater than or equal to that element
rows=ceil(sample_depth/el_size); 
mesh_size_x=sample_width/columns;
mesh_size_z=sample_depth/rows;
rows =rows+1;
columns =columns+1;

clear nodes
nodes(:,1)=repmat([0:mesh_size_x:sample_width],1,rows)'; % returns an array with n copies (repmat)
temp_z=ones(columns,1)*[1:rows];
temp_z=temp_z(:)';
nodes(:,2)=(temp_z-1)*mesh_size_z;
elements_temp = [1:columns-1;2:columns;columns+2:2*columns;columns+1:2*columns-1]';
elements = elements_temp;
for i=2:rows-1
    elements=cat(1,elements,elements_temp+columns*(i-1));
end

mod.nds = nodes;
mod.els = elements;

% when creating an element, be careful about how the elements are
% generated. 
rows =rows-1;
columns =columns-1;

mod.rows = rows;
mod.cols = columns;


%Now remove elements outside original boundary
% [in, out] = fn_elements_in_region(mod, bdry_pts);
% mod.els(out, :) = [];
% 
% %Tidy up by removing unused nodes
% [mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);

%Associate each element with a material index = 1
n_els = size(mod.els, 1);
mod.el_mat_i = ones(n_els, 1);
end
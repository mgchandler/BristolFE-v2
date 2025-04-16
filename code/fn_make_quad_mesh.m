function [nodes,elements,columns,rows] = fn_make_quad_mesh(element_size,sample_width,sample_depth)
tic

columns=ceil(sample_width/element_size); % ceil function rounds each element to nearest integer greater than or equal to that element
rows=ceil(sample_depth/element_size); 
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


% when creating an element, be careful about how the elements are
% generated. 
rows =rows-1;
columns =columns-1;
toc
end
% temp3
function [VX,VY,V,C]=fn_neper_grains(dataname,sample_width,sample_depth)

fid = fopen([dataname]);
data = textscan(fid,'%s','Delimiter','\n');
data = data{1};
% vertex
vertex_start = find(strcmp(data,'**vertex'));
vertex_num = str2num(data{vertex_start+1});
[~,vertex(:,1),vertex(:,2),~,~] = textread(dataname,'%f%f%f%f%f',vertex_num,...
    'headerlines', vertex_start+1,'delimiter',' ','emptyvalue',NaN);
% edge
edge_start=find(strcmp(data,'**edge'));
edge_num=str2num(data{edge_start+1});
[~,edge(:,1),edge(:,2),~] = textread(dataname,'%d%d%d%d',edge_num,...
    'headerlines', edge_start+1,'delimiter',' ','emptyvalue',NaN);
% face
face_start=find(strcmp(data,'**face'));
face_num=str2num(data{face_start+1});
face_cell=data(face_start+2:face_start+1+4*face_num);
for i=1:length(face_cell)
    face_irow=str2num(face_cell{i});
    face(i,1:length(face_irow))=face_irow;
end
face_vertices=face(([1:face_num]-1)*4+1,:);
face_vertices(:,[1,2])=[];
face_edge=abs(face(([1:face_num]-1)*4+2,:));
face_edge(:,[1])=[];
% 
VX(1,:)=vertex(edge(:,1),1)';
VX(2,:)=vertex(edge(:,2),1)';
VY(1,:)=vertex(edge(:,1),2)';
VY(2,:)=vertex(edge(:,2),2)';
V=vertex;
for i=1:size(face_vertices,1)
    face_vertices_irow=face_vertices(i,:);
    face_vertices_irow(face_vertices_irow==0)=[];
    C{i,1}=face_vertices_irow;
end
% 
WIDTH = max(max(VX));
DEPTH = max(max(VY));
scale1 = WIDTH/sample_width;
scale2 = DEPTH/sample_depth;
if scale1 ~= scale2
    error('There is somthing wrong with the settings of sample_width and sample_depth');
end
scale = scale1;
V = V/scale;
VX = VX/scale;
VY = VY/scale;

end



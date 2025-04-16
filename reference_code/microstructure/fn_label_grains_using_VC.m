function [grains]=fn_label_grains_using_VC(ele_cen,V,C)
tic
ele_cen_1=ele_cen(:,1);
ele_cen_2=ele_cen(:,2);
in_ind=1:length(ele_cen_1);
V_1=V(:,1);
V_2=V(:,2);
grains=cell(length(C),1);
for i=1:length(C)
%     fn_show_progress(i,length(C),'Assigning grains',1,1);
    if sum(V(C{i},:)==inf)>0
    elseif sum(V(C{i},:)==-inf)>0
    else
        in=inpolygon(ele_cen_1,ele_cen_2,V_1(C{i}),V_2(C{i}));% in=1 for elementcentre in poly made up with V(C{i})
        grains{i}=in_ind(in);
        if mod(i,1000)==0
            'grains{i}'
            i
            toc
            tic
        end
    end
end
end
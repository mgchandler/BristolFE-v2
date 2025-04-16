function CA = fn_euler_c(C0,a1,a,a2)

R11 = cos(a1).*cos(a2)-sin(a1).*sin(a2).*cos(a);
R21 = -cos(a1).*sin(a2)-sin(a1).*cos(a2).*cos(a);
R31 = sin(a1).*sin(a);
R12 = sin(a1).*cos(a2)+cos(a1).*sin(a2).*cos(a);
R22 = -sin(a1).*sin(a2)+cos(a1).*cos(a2).*cos(a);
R32 = -cos(a1).*sin(a);
R13 = sin(a2).*sin(a);
R23 = cos(a2).*sin(a);
R33 = cos(a);
%% 宏观晶粒弹性矩阵
for i=1:length(R11)
    % 原始文献里的方式，将原来坐标系下的弹性矩阵等效到现在坐标系下后的弹性矩阵
    R{i}=[R11(i)^2 R12(i)^2 R13(i)^2 2*R12(i)*R13(i) 2*R13(i)*R11(i) 2*R11(i)*R12(i)
        R21(i)^2 R22(i)^2 R23(i)^2 2*R22(i)*R23(i) 2*R23(i)*R21(i) 2*R21(i)*R22(i)
        R31(i)^2 R32(i)^2 R33(i)^2 2*R32(i)*R33(i) 2*R33(i)*R31(i) 2*R31(i)*R32(i)
        R21(i)*R31(i) R22(i)*R32(i) R23(i)*R33(i) R22(i)*R33(i)+R23(i)*R32(i) R21(i)*R33(i)+R23(i)*R31(i) R22(i)*R31(i)+R21(i)*R32(i)
        R31(i)*R11(i) R32(i)*R12(i) R33(i)*R13(i) R12(i)*R33(i)+R13(i)*R32(i) R13(i)*R31(i)+R11(i)*R33(i) R11(i)*R32(i)+R12(i)*R31(i)
        R11(i)*R21(i) R12(i)*R22(i) R13(i)*R23(i) R12(i)*R23(i)+R13(i)*R22(i) R13(i)*R21(i)+R11(i)*R23(i) R11(i)*R22(i)+R12(i)*R21(i)];
    CA{i}=R{i}*C0*R{i}';
end


end
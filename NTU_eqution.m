%% Design of NTU eqution
%NTU=a*V*Fa^b*Fz^c*(1-ksi)^d
% function NTU=NTU_eqution(V,Fa,Fz,ksi,a,b,c,d)
function NTU=NTU_eqution(a,X)
% tc=X(:,1);
% te=X(:,2);
% C=a(1)*te+a(2)*te.^2+a(3)*te.^3+a(4)*te.*tc+a(5)*te.^2.*tc+a(6)*te.*tc.^2+a(7)*tc+a(8)*tc.^2+a(9)*tc.^3+a(10);
V=X(:,1);Fa=X(:,2);Fz=X(:,3);ksi=X(:,4);
NTU=a(1)*V*Fa^a(2)*Fz^a(3)*(1-ksi)^a(4);
end
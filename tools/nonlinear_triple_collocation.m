function [NE,NI] = nonlinear_triple_collocation(X,Y,Z,B)

% basic stats
[Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,B,B,B);

% correlations
NI(1) = (Ixy+Ixz-Ixyz)/Hx;   
NI(2) = (Ixy+Iyz-Ixyz)/Hy;   
NI(3) = (Ixz+Iyz-Ixyz)/Hz;

% residuals
NE(1) = 1-(Ixy+Ixz-Ixyz)/Hx; 
NE(2) = 1-(Ixy+Iyz-Ixyz)/Hy; 
NE(3) = 1-(Ixz+Iyz-Ixyz)/Hz;



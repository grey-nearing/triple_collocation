function [sX,sY,sZ,rX,rY,rZ] = triple_collocation(X,Y,Z)

Qxx = cov(X,X); Qxx = Qxx(2);
Qyy = cov(Y,Y); Qyy = Qyy(2);
Qzz = cov(Z,Z); Qzz = Qzz(2);

Qxy = cov(X,Y); Qxy = Qxy(2);
Qxz = cov(X,Z); Qxz = Qxz(2);
Qyz = cov(Y,Z); Qyz = Qyz(2);

sX = Qxx - Qxy*Qxz/Qyz;
sY = Qyy - Qxy*Qyz/Qxz;
sZ = Qzz - Qxz*Qyz/Qxy;

sX = sX/Qxx;
sY = sY/Qyy;
sZ = sZ/Qzz;

rX = Qxy*Qxz/Qxx/Qyz;
rY = Qxy*Qyz/Qyy/Qxz;
rZ = Qxz*Qyz/Qzz/Qxy;



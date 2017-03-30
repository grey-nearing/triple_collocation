function [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,Bx,By,Bz)

N = length(X);
assert(N==length(Y))
assert(N==length(Z))

Mx = length(Bx);
My = length(By);
Mz = length(Bz);

Pxyz = zeros(Mx,My,Mz);
for n = 1:N
  x = X(n);
  y = Y(n);
  z = Z(n);
  y = find(y<By,1,'first');
  x = find(x<Bx,1,'first');
  z = find(z<Bz,1,'first');
  if isempty(x); fprintf('X: min=%f, max=%f, Bmin=%f, Bmax=%f\n',min(X),max(X),Bx(1),Bx(end)); keyboard; end
  if isempty(y); fprintf('Y: min=%f, max=%f, Bmin=%f, Bmax=%f\n',min(Y),max(Y),By(1),By(end)); keyboard; end
  if isempty(z); fprintf('Z: min=%f, max=%f, Bmin=%f, Bmax=%f\n',min(Z),max(Z),Bz(1),Bz(end)); keyboard; end
  Pxyz(x,y,z) = Pxyz(x,y,z) + 1;
end
Pxyz = Pxyz/N;

Pxy = squeeze(sum(Pxyz,3));
Pxz = squeeze(sum(Pxyz,2));
Pyz = squeeze(sum(Pxyz,1));
Px  = squeeze(sum(Pxy,2));
Py  = squeeze(sum(Pxy,1));
Pz  = squeeze(sum(Pxz,1));

Px   = Px(:);
Py   = Py(:);
Pz   = Pz(:);
Pyz  = Pyz(:);
Pxz  = Pxz(:);
Pxy  = Pxy(:);
Pxyz = Pxyz(:);

if abs(sum(Px)-1)>1/N^2; error('Px does not sum to 1'); end;
if abs(sum(Py)-1)>1/N^2; error('Py does not sum to 1'); end;
if abs(sum(Pz)-1)>1/N^2; error('Pz does not sum to 1'); end;

if abs(sum(Pxy)-1)>1/N^2; error('Pxy does not sum to 1'); end;
if abs(sum(Pxz)-1)>1/N^2; error('Pxz does not sum to 1'); end;
if abs(sum(Pyz)-1)>1/N^2; error('Pyz does not sum to 1'); end;

if abs(sum(Pxyz)-1)>1/N^2; error('Pxyz does not sum to 1'); end;

Hx = -Px(Px>0)'*log(Px(Px>0));
Hy = -Py(Py>0)'*log(Py(Py>0));
Hz = -Pz(Pz>0)'*log(Pz(Pz>0));

Hxy = -Pxy(Pxy>0)'*log(Pxy(Pxy>0));
Hxz = -Pxz(Pxz>0)'*log(Pxz(Pxz>0));
Hyz = -Pyz(Pyz>0)'*log(Pyz(Pyz>0));

Hxyz = -Pxyz(Pxyz>0)'*log(Pxyz(Pxyz>0));

%I = Hxz+Hxy-Hxyz-Hx;
Ixy = Hx+Hy-Hxy;
Ixz = Hx+Hz-Hxz;
Iyz = Hy+Hz-Hyz;
Ixyz = -(Hx+Hy+Hz-Hxyz-Ixy-Ixz-Iyz);

Ixy_z = Ixz + Hxz + Hyz - Hxyz - Hz;




function [Ixy,Hx,Hy] = mutual_info(X,Y,Bx,By)

N = length(X);
assert(N==length(Y))

Mx = length(Bx);
My = length(By);

Pxy = zeros(Mx,My);
for n = 1:N
  x = X(n);
  y = Y(n);
  y = find(y<By,1,'first');
  x = find(x<Bx,1,'first');
  if isempty(x) || isempty(y); display(n); keyboard; end
  Pxy(x,y) = Pxy(x,y) + 1;
end
Pxy = Pxy/N;

Px  = squeeze(sum(Pxy,2));
Py  = squeeze(sum(Pxy,1));

Pxy  = Pxy(:);
Px   = Px(:);
Py   = Py(:);

if abs(sum(Px)-1)>1/N^2; error('Px does not sum to 1'); end;
if abs(sum(Py)-1)>1/N^2; error('Py does not sum to 1'); end;
if abs(sum(Pxy)-1)>1/N^2; error('Pxy does not sum to 1'); end;

Hx = -Px(Px>0)'*log(Px(Px>0));
Hy = -Py(Py>0)'*log(Py(Py>0));

Hxy = -Pxy(Pxy>0)'*log(Pxy(Pxy>0));

Ixy = Hx+Hy-Hxy;




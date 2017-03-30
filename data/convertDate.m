function doy = convertDate(date)

y = date(1);
m = date(2);
d = date(3);
h = date(4);
mn = date(5);

N = [31,28,31,30,31,30,31,31,30,31,30,31];
if mod(y,4) == 0; N(2) = 29; end;

doy = 0;
if m > 1
 doy = sum(N(1:m-1));
end

doy = doy + d + h/24 + mn/24/60;



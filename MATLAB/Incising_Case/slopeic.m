function u0 = slopeic(x)
a = 0;
b = 0.02;
c = 0.01;
d = 0;
e = 0;
u0 = a.*x.^(3) + b.*x.^2 + c.*x + d;
end
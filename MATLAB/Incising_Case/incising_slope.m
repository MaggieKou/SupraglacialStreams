function [c,f,s] = incising_slope(x, t, u, dudx)
c = 1;
f = u.^(6/5)./((1+u.^2).^(8/5));
s = 0;
end
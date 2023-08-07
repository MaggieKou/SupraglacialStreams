function y = sinusoidal(x, a, b, c, d, e)
y = real(a*sin(b.*(x-c)) + d + e.*x);
end
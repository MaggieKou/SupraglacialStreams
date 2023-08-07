function y = tanh_fun(x, a, b, c, d, e)
y = real(a*tanh(b.*(x-c)) + d + e.*x);
end
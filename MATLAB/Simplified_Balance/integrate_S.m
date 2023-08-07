function h = integrate_S (h_start, S_vec, i, dx)
h = h_start;
for k = 1:i
    h = h + S_vec(k)*dx; %riemann sum over previous values 
end
end
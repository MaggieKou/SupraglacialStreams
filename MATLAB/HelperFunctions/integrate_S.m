function h = integrate_S (dx, Sj, k, h_o)
h = h_o;
for i = 1:k
    h = h + dx*Sj(i); %integrator 
end
end
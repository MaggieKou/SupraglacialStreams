function w = integrator(dx, u, x, k, h_o)
%u is given per time-step, is also the same length as x
w = h_o;
for i = 1:k
    w = w + dx/sqrt((u(i)^(10) -1)); %integrator 
end
end

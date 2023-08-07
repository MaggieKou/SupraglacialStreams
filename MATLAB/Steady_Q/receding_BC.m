function [h_x_0, S_0] = receding_BC(h_init, S_init,nt, v, dt, Q)
%given scalars h_in and S_in, we march forward as an ODE with S_t = 5/2vS(1-S^2)^4/5(sqrt(1+S^2) + 1/v*S)
h_x_0 = zeros([1, nt]);
h_x_0(1) = h_init;
S_0 =  zeros([1, nt]);
S_0(1) = S_init;

for j = 2:nt
    S_0 (j) = S_0(j-1) + 5/2*Q^(-2/5)*dt*(S_0(j-1)^(11/5)*(1+S_0(j-1)^2)^(3/10)*((1 + S_0(j-1)^2)^(3/5)*v - S_0(j-1)^(1/5)*Q^(3/5)));
    h_x_0(j) = h_x_0(j-1) + dt*v*S_0(j);
end
end
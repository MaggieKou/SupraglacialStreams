function [H_pp,out] = receding_solver(S, v, dt)
%given S, input slope at a given time step, j
H_p = v^5.*S;

H_pp = 5/2.*(1+v^(10).*H_p.^2).^(3/10).*H_p.^(11/5).*(-H_p.^(1/5)+(1+v^(10).*H_p.^2).^(3/5));
out = S + dt*v^(12).*H_pp;
end
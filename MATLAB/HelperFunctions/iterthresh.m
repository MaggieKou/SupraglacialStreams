function [uret,h, num] = iterthresh(err, v, hj, Qj, Qtj, uj, dx, x, dt, h_x_0)
num = 0; %function counts number of iterations until threshold is reached
uret = uj;
h = zeros(length(uret));
u_new  = init_guess(v, dt, Qj, Qtj, uj);

%find new value on interation
for i = 1:length(x)
%initial guess for the integral at each point
w_new = integrator(dx, uret, x, i, h_x_0(i));
der_int = derivative_calc(w_new, hj(i), dt);
unew2 = pde_solve(Qj(i), Qtj(i), der_int, uj(i), dt);
%iterate until we reach the error threshold
while (abs(max(unew2 - u_new)) > err) 
u_new = unew2;
w_new = integrator(dx, uret, x, i, h_x_0(i)); %calculates h at the current timestep

der_int = derivative_calc(w_new, hj(i), dt); %calcs derivative with new value 
unew2 = pde_solve(Qj(i), Qtj(i), der_int, uj(i), dt); %keeps solving with new value
num = num +1;
end
h(i) = w_new;
uret(i) = unew2;
end
end
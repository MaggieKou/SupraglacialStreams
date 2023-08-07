function [Sj2, h2] = S_solver(err, Qj, Qtj, Sj, v, dt,  dx, h_x_0)
%given the slope at a previous timestep, returns slope at later timestep

h2 = zeros(length(Sj));
%initial guess
Sj_guess = Sj + dt.*(5/2.*Qj.^(-1/5).*(1+Sj.^2).^(3/10).*Sj.*(4/5.*Qtj.*(1+Sj.^2).^(7/10) - Sj.^(7/5).*Qj.^(3/5)));
%find new value on interation
for i = 1:length(Qj)
%initial guess for the integral at each point
h_cur = integrate_S(dx, Sj, i, h_x_0(i));
h_new = integrate_S(dx, Sj_guess, i, h_x_0(i));

der_int = derivative_calc(h_new, h_cur, dt); %calculates derivative 

Sj2 = Sj + dt.*(5/2.*Qj.^(-1/5).*(1+Sj.^2).^(3/10).*Sj.*(4/5.*Qtj.*(1+Sj.^2).^(7/10)+Sj.^(1/5).*(1+Sj.^2).^(3/5).*Qj.^(3/5).*(der_int) - Sj.^(7/5).*Qj.^(3/5)));
%iterate until we reach the error threshold
while (abs(max(Sj2 - Sj_guess)) > err) 
    %update from last iter
Sj_guess = Sj2;
h_cur = h_new; 

%new iter
h_new = integrate_S(dx, Sj_guess, x, i, h_x_0(i));
der_int = derivative_calc(h_new, h_cur, dt);
%update
Sj2 = Sj + dt.*(5/2.*Qj.^(-1/5).*(1+Sj.^2).^(3/10).*Sj.*(4/5.*Qtj.*(1+Sj.^2).^(7/10)+Sj.^(1/5).*(1+Sj.^2).^(3/5).*Qj.^(3/5).*(der_int) - Sj.^(7/5).*Qj.^(3/5)));
end
h2(i) = h_new;
end
end
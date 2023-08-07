function u_new = pde_solve(Qj, Qtj, der_int, uj, dt)
    u_new=  uj + dt/2.*(Qj.^(1/5).*(uj).^(-6) - 4/5.*(Qj).^(-1).*Qtj.*uj - Qj.^(-2/5).*der_int);
end

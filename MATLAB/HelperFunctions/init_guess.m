function u_new = init_guess(v, dt, Qj, Qtj, uj)
%assume the different is zero for integrating
ut = (Qj.^(1/5).*(uj).^(-6)- 4/5.*(Qj).^(-1).*Qtj.*uj- Qj.^(-2/5).*((v+Qtj).*(uj.^(10) -1).^(-1/2)));
u_new = uj + dt/2.*ut;
%- Qj.^(-2/5).*(uj.^(10) -1).^(-2)
end
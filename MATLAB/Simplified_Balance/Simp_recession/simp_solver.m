function h_eta = simp_solver(Q, h_eta_init, t)
h_eta = zeros([1,length(t)]);
dt = (t(length(t))-t(1))/length(t);
h_eta(1) = h_eta_init;

for i = 2:length(t)
    h_eta(i) = h_eta(i-1) + dt*(1e-7)*(5/2.*Q.^(-2/5)*h_eta(i-1).*(1- h_eta(i-1).^(6/5))- h_eta(i-1).^(-1));
end
end
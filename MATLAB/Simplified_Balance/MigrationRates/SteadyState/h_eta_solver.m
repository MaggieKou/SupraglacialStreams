function h_eta_sol = h_eta_solver(n_e, n_t, d_t, a, d_eta, t, Q, hx_init)
%prescribe bcs
h_eta_0 = max(hx_init); %slope threshold from input data

h_eta_sol = zeros([n_t, n_e]);
h_eta_sol(:,1) = h_eta_0; %BC from max slope in travelling reference frame
h_eta_sol(1, :) = hx_init; %IC from initial slope

for j = 2:n_t
    for i = 1: n_e
        h_eta_sol(j,i) =  h_eta_sol(j-1, i) + 2/a*t(j)^(-1/2)*d_t*5*h_eta_sol(j-1, i)^(11/5)*Q^(-2/5)*a*t(j)^(-1/2)*(h_eta_sol(j-1, i)^(1/5)*Q^(3/5)+1/2*t(j)^(-1/2)*a);
    end
end

h_eta_2 = zeros([n_t, n_e]);

for j = 2:n_e
    h_eta_2(:, j) = (h_eta_sol(:, j) - h_eta_sol(:,j-1))/d_eta; 
end 
end

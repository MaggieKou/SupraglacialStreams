addpath 'C:\Users\maggi\Documents\School Stuff\ESROP\Code\SupraglacialStreams\MATLAB\HelperFunctions'

n_e = 1000; 
eta_o = 0;
eta_f = 0.8;
eta = linspace(eta_o, eta_f, n_e);
d_eta = (eta_f-eta_o)/n_e;

n_t = 10000; 
t_o = 1e-1; 
t_f = 1700; 
t = linspace(t_o, t_f, n_t);
d_t = (t_f-t_o)/n_t;

Q = 1e-4;
a_fit = 0.290; %from fit

%IC: 
a = 1*10^(-3);
b = 8;
c = 0.5;
d = 0.01; 
e = 1.*10^(-10);

h_init = tanh_fun(eta, a, b, c, d, e);
hx_init = zeros([1, n_e]);
for j = 2:n_e
    hx_init(j) = (h_init(j) - h_init(j-1))/d_eta; 
end
hx_init(1) = hx_init(2);

h_eta_out = h_eta_solver(n_e, n_t, d_t, a_fit, d_eta, t, Q, hx_init); 
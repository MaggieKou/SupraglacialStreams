addpath 'C:\Users\maggi\Documents\School Stuff\ESROP\Code\SupraglacialStreams\MATLAB\HelperFunctions'
%initialization 
nx = 500;
xo = 0;
xf = 1;
dx = (xf-xo)/nx;

nt = 1700; 
tf = 1700;
dt = tf/nt;

x = linspace(xo,xf,nx);
t = linspace(0,tf,nt);

% %prescribe initial profile 
% a = 0.000000;
% b = 1.*10^(-5);
% c = 5.*10^(-5);
% d = 0; 

%tanh profile
a = 1e-2;
b = 8;
c = 0.5;
d = 1e-2;
e = 1e-2;

%sinusoidal profile
% a = 0.000004;
% b = 8;
% c = 0.5;
% d = 0.01; 
% e = 0.000005;

Q = 1*10^(-9);

h = zeros([nt, nx]);
h_init = tanh_fun(x, a, b, c, d, e); %y = real(a*tanh(b.*(x-c)) + d + e.*x);
%h_init = poly_fun(x, a, b, c, d);
%h_init = sinusoidal(x, a, b, c, d, e);

h(1, :) = h_init;
S = zeros([nt, nx]);
%initial slope
for i = 2:nx
    S(1, i) = (h(1, i) - h(1, i-1))/dx;
end
S(1,1) = S(1, 2); %assume smoothness


for j = 1: nx
    S(:, j) = simp_solver(Q, S(1, j), t);
end

for j = 2:nt
    h(j, 1) = h(j-1, 1) + dt*S(j-1, 1).^(22/5)*Q^(3/5)/(2/5*Q^(2/5) + S(j-1, 1)^(16/5));
end 

for i = 2: nx
    h(:, i) = h(:, i-1) + dx.*S(:, i-1);
end

h_xt =zeros([nt, nx]);
%calculate h_xt
for j = 2:nt
    h_xt(j,:) = (S(j, :) - S(j-1, :))/dt;
end
h_xt(1,:) = h_xt(2,:);

v = (S.^(17/5)*Q^(3/5))./(2/5*Q^(2/5) + S.^(16/5));
eta = (x + t*v)./(2*x - 1);

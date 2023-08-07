addpath 'C:\Users\maggi\Documents\School Stuff\ESROP\Code\SupraglacialStreams\MATLAB\HelperFunctions'
%initialization 
nx = 500;
xo = 0;
xf = 1;
dx = (xf-xo)/nx;

nt = 1000; 
tf = 1700;
dt = tf/nt;

x = linspace(xo,xf,nx);
t = linspace(0,tf,nt);

% %prescribe initial profile 
a = 0.000000;
b = 0.00000;
c = 0.00001;
d = 0; 
e = 0.01;

%tanh profile
% a = 0.000004;
% b = 8;
% c = 0.5;
% d = 0.01; 
% e = 0.000005;

%sinusidal profile
% a = 0.000004;
% b = 8;
% c = 0.5;
% d = 0.01; 
% e = 0.000005;

v = 4*10^(-10); 
Q = 8*10^(-7);
err = 0.001;

%h_init = tanh_fun(x, a, b, c, d, e); %y = real(a*tanh(b.*(x-c)) + d + e.*x);
h_init = poly_fun(x, a, b, c, d);
%h_init = sinusoidal(x, a, b, c, d, e);
h = zeros([nt, nx]);

% x_shift = x;
% for i = 1:nt
%     x_shift(i) = x(1)+ v*t(i);
% end
% h_x_0 = poly_fun(x_shift, a, b, c, d);
% %h_x_0 = tanh_fun(x_shift, a, b, c, d, e); %limits definition for h_init
% %h_x_0 = sinusoidal(x_shift, a, b, c, d, e);
h(1, :) = h_init;
S = zeros([nt, nx]);
%initial slope
for i = 2:nx
    S(1, i) = (h(1, i) - h(1, i-1))/dx;
end
S(1,1) = S(1, 2); %assume smoothness
[h_x_0, S(:,1)]= receding_BC(h_init(1),S(1,1),nt, v, dt, Q);

% %time stepping scheme
for j = 2:nt
    h(j, :) = row_iterate(h_x_0(j), S(j-1, :), dx, nx, dt, Q, err);
    for i = 2:nx
    S(j, i) = (h(j, i)-h(j, i-1))/dx;
    end
    S(j, 1) = S(j, 2);
end

Sx=S; 
St = zeros([nt,nx]);
ht = zeros([nt,nx]);
Sxt = zeros([nt,nx]);
for i = 2:nx
Sx(:,i) = (S(:, i)-S(:, i-1))./dx;
end
Sx(:, 1) = S(:, 2);
Sx(2, :) = Sx(3, :);

for j = 2:nt
ht(j, :) = (h(j,:)-h(j-1,:))./dt;
St(j, :) =(S(j, :)-S(j-1, :))./dt;
Sxt(j, :)=(Sx(j, :)-Sx(j-1, :))./dt;
end
St(1, :) = St(2, :);
ht(1, :) = ht(2, :);
Sxt(1, :) = Sxt(2, :);
Sxt(:, 2) = Sxt(:, 3);

w = S.^(-1/5);
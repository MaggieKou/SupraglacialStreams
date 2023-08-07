addpath 'C:\Users\maggi\Documents\School Stuff\ESROP\Code\SupraglacialStreams\MATLAB\HelperFunctions'

%initialization 
nx = 50;
xo = 0;
xf = 10;
dx = (xf-xo)/nx;

nt = 1000; 
tf = 17000;
dt = tf/nt;

x = linspace(xo,xf,nx);
t = linspace(0,tf,nt);

%prescribe initial profile 
a = 0.00;
b = 0.0001;
c = 0.00001;
d = 0; 
e = 0.01;

Q = 10^(-6);
%h_init = tanh_fun(x, a, b, c, d, e); %y = real(a*tanh(b.*(x-c)) + d + e.*x);
h_init = poly_fun(x, a, b, c, d);
h = zeros([nt, nx]);
h(1, :) = h_init;
S = zeros([nt, nx]);
%initial slope
for i = 2:nx
    S(1, i) = (h(1, i) - h(1, i-1))/dx;
end
S(1,1) = S(1, 2); %assume smoothness

% %time stepping scheme
for j = 2:nt
    h(j, 1) = h(j-1, 1) + dt*Q^(3/5)*(S(j-1,1)^(6/5))/(1+S(j-1,1)^2)^(3/5); %dynamically updates BC for receding case
    for i = 2: nx
    h(j, i) = h(j-1, i) + dt*Q^(3/5)*(h(j-1, i)-h(j-1, i-1))^(6/5)/((dx^2 + (h(j-1, i)-h(j-1, i-1))^2)^(3/5)); %dynamically updates
    S(j, i) = (h(j,i) - h(j, i-1))/dx;
    end
    S(j, 1) = S(j, 2); 
end

Sx=S; 
St = zeros([nt,nx]);
for i = 2:nx
Sx(:,i) = (S(:, i)-S(:, i-1))./dx;
end

for j = 2:nt
St(j, :) =(S(j, :)-S(j-1, :))./dt;
end

% %space step upwinding scheme 
% for j = 2:nt
%     h(j, 1) = h(j-1, 1) + dt*Q^(3/5)*(S(j-1,1)^(6/5))/(1+S(j-1,1)^2)^(3/5); %dynamically updates BC for receding case
%     for i = 2: nx
%     h(j, i) = h(j, i-1) + dx/(((h(j,i-1)-h(j-1, i-1))/dt)^(-5/3)*Q-1)^(1/2); %dynamically updates
%     S(j, i) = (h(j,i) - h(j, i-1))/dx;
%     end
%     S(j, 1) = S(j, 2); 
% end

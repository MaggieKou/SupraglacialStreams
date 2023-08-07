%given our initial condition S(t = 0, x), we can derive u = ((1 +
%S^2)^1/2/S)^1/5 at t = 0 and integrate%

%SET UP THE SOLUTION GRID
nx = 100;
xo= 0;
xf= 1;
x = linspace(xo, xf, nx);
dx = (xf-xo)/nx;

nt = 17000;
to = 0;
tf = 1700;
t = linspace(to, tf, nt);
dt = (tf-to)/nt;

%error threshold 
err = 0.001;

%need to define Q(t, x) and u_initial

%import data from Q timeseries
T = readtable("RawData/Discharge_timeseries.csv");
Q = zeros([nt, nx]);
tau = 27.3;
for i = 1:nx
    %Q(:, i) = 10^(-4).*table2array(T(1:17000,5));
    for j = 1:nt
    %Q(j, i) = 10^(-9).*(1 - 0.9*sin(2*pi/tau.*t(j)));
    Q(j, i) =10^(-9);
    end
end
Qt = zeros([nt, nx]);
for i = 1:nx
    for j = 2:nt
        Qt(j,i) = (Q(j,i) - Q(j-1,i))/dt;
    end
end
Qt(1,:) = Qt(2,:);

u  = zeros([nt,nx]);
a = 0.00;
b = 0.1; %0.0005 yielded convergence to steady state slope 
c = 0.01;
d = 0;
e = 0.0001;

v =10^(-5);  %recession velocity, observed final convergent slopes around 10^-10 to 10^-8 ranges around an order of mag less than observed

h_init =poly_fun(x, a, b, c, d);
%h_init = tanh_fun(x, a, b, 
% c, d, e);

x_shift = x;
for i = 1:nt
    x_shift(i) = x(1)- v*t(i);
end

%h_x_0 = poly_fun(x_shift, a, b, c, d);
%h_x_0 = tanh_fun(x_shift, a, b, c, d, e); %limits definition for h_init;
h = zeros([nt,nx]);
h2 = zeros([nt,nx]);
S = zeros([nt,nx]);

h(1,:) = h_init;
h2(1, :) = h_init;

for i = 2:nx
S(1, i) = (h_init(i)- h_init(i-1))/dx;
u(1,i) = ((1+S(1,i)^2)^(1/2)/S(1,i))^(1/5);
end
S(1,1) = S(1,2);
h_rec = receding_method(h, nx, dx, nt, dt, v);
u(1,1) = u(1, 2); 
%go through solution, step by step
for j = 2:nt
    u(j,:) = iterthresh(err,v, h(j-1, :), Q(j,:), Qt(j,:), u(j-1,:),dx, x, dt, h_rec);
    %u(j,:) = init_guess(v, dt, Q(j,:), Qt(j,:), u(j-1,:));
    %S(j,:); h2(j, :) = S_solver(err, Q(j, :), Qt(j,:), S(j-1, :), v, dt, dx, h_rec);
end

%now back-calculate slope 
for j = 1:nt
    for i = 1:nx
        S (j,i) = 1/(sqrt(u(j,i)^(10)-1));
    end
end

%h(:, 1) = h_x_0;

for i = 2:nx
    h(:,i) = h(:,i-1) + dx.*S(:,i);
    %h(j, :) = h(j-1,:) + dt.*(u(j,:).^(-6).*Q(j,:).^(3/5) + 2*(Q(j,:).^(2/5).*u(j,:) - Q(j-1,:).^(2/5).*u(j-1,:))/dt);
end


w = ((1+S.^2)./S).^(1/5); %width
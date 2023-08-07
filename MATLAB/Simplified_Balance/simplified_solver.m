addpath 'C:\Users\maggi\Documents\School Stuff\ESROP\Code\SupraglacialStreams\MATLAB\HelperFunctions'
%initialization 
nx = 1000;
xo = 0;
xf = 0.8;
dx = (xf-xo)/nx;

to = 100;
nt = 100; 
tf = 200; %dt = 17 is like half a day resolution for each timestep
dt = tf/nt;

x = linspace(xo,xf,nx);
t = linspace(to,tf,nt);

% %prescribe initial profile 
% a = 0.000000;
% b = 1.*10^(-5);
% c = 5.*10^(-5);
% d = 0; 

%tanh profile
a = 5e-3;
b = 8;
c = 0.5;
d = 5e-3;
e = 5e-3;

%sinusoidal profile
% a = 0.000004;
% b = 8;
% c = 0.5;
% d = 0.01; 
% e = 0.000005;

v = 10^(-6); 
Q = 10^(-2);

err = 0.001;

h_init = tanh_fun(x, a, b, c, d, e); %y = real(a*tanh(b.*(x-c)) + d + e.*x);
%h_init = poly_fun(x, a, b, c, d);
%h_init = sinusoidal(x, a, b, c, d, e);
h = zeros([nt, nx]);
h_x_0 = zeros([nt,1]);
h_x_0(1) = h_init(1);

%h_x_0 = poly_fun(x_shift, a, b, c, d);
%h_x_0 = tanh_fun(x_shift, a, b, c, d, e); %limits definition for h_init
%h_x_0 = sinusoidal(x_shift, a, b, c, d, e);

h(1, :) = h_init;
S = zeros([nt, nx]);
%initial slope
for i = 2:nx
    S(1, i) = (h(1, i) - h(1, i-1))/dx;
end
S(1,1) = S(1, 2); %assume smoothness
[h_x_0, S(:,1)]= receding_BC(h_init(1),S(1,1),nt, v, dt, Q);


%solves for the boundary condition on the assumption that the slope would
%be kept the same for all t
for i = 2:nt
    h_x_0(i) = h_x_0(i-1)+ S(1,1)^(6/5)*Q^(3/5)*dt;
end

% %time stepping scheme
for j = 2:nt
    h(j, :) = row_iterate(h_x_0(j), h_x_0(j-1), S(j-1, :), dx, nx, dt, Q, err); 
    for i = 2:nx
    S(j, i) = (h(j, i)-h(j, i-1))/dx;
    end
    S(j, 1) = S(j, 2);
end

Sx=S; 
St = zeros([nt,nx]);
ht = zeros([nt,nx]);
htt = zeros([nt,nx]);
Sxt = zeros([nt,nx]);
%calculate incision rate and values
for i = 2:nx
Sx(:,i) = (S(:, i)-S(:, i-1))./dx;
end
Sx(:, 1) = S(:, 2);
Sx(2, :) = Sx(3, :);

for j = 2:nt
ht(j, :) = (h(j,:)-h(j-1,:))./dt;
htt(j,:) =  (ht(j,:)-ht(j-1,:))./dt;
St(j, :) =(S(j, :)-S(j-1, :))./dt;
Sxt(j, :)=(Sx(j, :)-Sx(j-1, :))./dt;
end
St(1, :) = St(2, :);
ht(1, :) = ht(2, :);
htt(1, :) = htt(2, :);
Sxt(1, :) = Sxt(2, :);
Sxt(:, 2) = Sxt(:, 3);

prop_speed = Sx./(1/5.*S.^(-6/5).*Q^(2/5)+sqrt(1/25.*S.^(-12/5).*Q.^(-4/5)+Sx.*htt));
prop_speed2 = 5/6.*S.^(1/5).*Q^(3/5);
slope_thresh = 0.98*max(S(1,:));
RHm = zeros([1, nt]); %rhs migration rate
LHm = zeros([1, nt]); %lhs migration rate
v_recR = zeros([1, nt]);
v_recL = zeros([1, nt]);
St_l = zeros([1, nt]);
ht_l = zeros([1, nt]);
St_r = zeros([1, nt]);
ht_r = zeros([1, nt]);

%find max in the propagation speed: 
x_vp_max = zeros([1, nt]);
x_vp2_max = zeros([1, nt]);
for j = 1:nt
    vp_max_temp = max(prop_speed(j,:));
    vp2_max_temp = max(prop_speed2(j,:));
    for i = 1:nx
        if prop_speed(j, i) == vp_max_temp
            x_vp_max(j) = x(i);
        end 
        if prop_speed2(j,i) == vp2_max_temp
            x_vp2_max(j) = x(i);
        end
    end
end

for i = 1:5
[RHm(i,:), LHm(i,:), v_recR(i,:), v_recL(i,:), St_r(i,:), St_l(i,:), ht_r(i,:), ht_l(i,:)] = reced_tracker(x, S, St, ht, (1-0.005*i)*max(S(1,:)), nt, nx, dt, dx); %call function to track slope, etc...
%single max thresh [RHm, LHm, v_recR, v_recL, St_r, St_l, ht_r, ht_l] = reced_tracker(x, S, St, ht, slope_thresh, nt, nx, dt, dx); %call function to track slope, etc...
end

%plotting in moving reference
eta = zeros(5, nt, nx);

for i = 1:5
for j = 1:nt
eta(i, j, :) = x - LHm(i, j);
end
end

eta_p = zeros([5,nt]);
S_eta_p = zeros([5,nt]);

eta_norm = zeros(5, nt, nx);

for i = 1:5
for j = 1:nt
    eta_norm(i,j,:) = eta(i,j,:)./(-RHm(i, j)+ LHm(i, j));
end
end

%new portion of the code takes forever to run, please try to time/ make
%more efficient
for i = 1:5
    for j = 1:nt
    for k = 1:nx
        if S(j,k) == (1-0.005*i)*max(S(j,:))
        S_eta_p(i,j) = S(j,k);
        eta_p(i,j) = eta(i,j,k);
        else if k>1 && (S(j, k-1) - (1-0.005*i)*max(S(j,:)))*(S(j,k)-(1-0.005*i)*max(S(j,:))) <=0
        S_eta_p(i,j) = S(j,k);
        eta_p(i,j) = eta(i,j,k);
        end
        end
    end
    end
end



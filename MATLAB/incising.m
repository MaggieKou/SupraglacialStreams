%initialization 
x = linspace(0,10,1000);
t = linspace(0,10,1000);

[sol,tsol,sole,te,ie] = pdepe(0,@incising_sol,@riveric,@riverbc,x,t);

function [c,f,s] = incising_sol(x,t,u,dudx)
c = -1;
f = 0;
s = (dudx)^(6/5)/(sqrt(1+(dudx)^2));
end

function h0 = riveric(x) % initializes the initial conditions for h
a = 0;
b = -2;
c = 0;
d = 50;
h0 = a*x^3 + b*x^2 + c*x + d;
end

function [pl, ql, pr, qr] = riverbc(xl,ul,xr,ur,t) %initializes the boundary conditions for h
Cf = 0.1;
g = 9.8;
Q = 4;
R = 0.5;
%BC
wt = linspace(1, 0, 1000); %rate of width changing with time is linearly decreasing

w0 = -(g*(ur-ul)/((xr-xl)*Cf*R*Q^2))^(-1/5); %initial width per BC
w = zeros(length(wt));
w(1) = w0;

%initializes w
for i = 2:length(wt)
    w(i) = w(i-1) + ((t(length(t)) -t(1))/length(t))*wt(i);
end

S_bc = Cf*R*Q^2./(g*w.^5); 
h_bc = zeros(length(wt));
h_bc(1) = ul(1);
for i = 2:length(h_bc)
    h_bc(i) = h_bc(i-1) + ((t(length(t)) -t(1))/length(t))*(S_bc(i+1)+S_bc(i-1))/2;
end
ql = 0;
qr = 0;
pr = ur - h_bc(length(t)); %sets boundary conditions along the domain
pl = ul - h_bc(1);
end

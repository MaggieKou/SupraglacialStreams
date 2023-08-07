
function ret = receding_method(h,nx, dx, nt, dt,v)

ret = zeros([nt, nx]);
ret(1, :) = h(1,:); % IC
S = zeros([nt, nx]);
H_pp = zeros([nt, nx]);

for i = 2:nx
        S(1, i) = (h(1, i) - h(1, i-1))/dx;
end
S(1,1) = S(1, 2);

for j = 2:nt
    ret(j, :) = ret(j-1, :) + dt*v*S(j-1,:); %takes dh/dt as a function of the slope and recession velocity 
    [H_pp(j,:), S(j, :)] = receding_solver(S(j-1, :), v, dt);
end

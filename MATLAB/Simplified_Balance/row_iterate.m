function h = row_iterate(h_start, h_b4, S_vec, dx, nx, dt, Q, err)
%start with previous row 
h_diff = (h_start - h_b4)./dt;
S = pde_Solve(S_vec, Q, dt, h_diff);
h = zeros([1, length(S)]);
dhdt = zeros([1,length(S)]);
%start with integrating the slope along the line 
for i = 1:nx
    dhdt(i) = (integrate_S(h_start, S_vec, i, dx) - h(i))/dt;
    h(i) = integrate_S(h_start, S_vec, i, dx); % new h value 
end

%now recalcuate with pde
S_guess = pde_Solve(S, Q, dt, dhdt);

while (max(abs(S_guess - S)) > err)
    S = S_guess;
    for i = 1:nx
    dhdt(i) = (integrate_S(h_start, S_vec, i, dx) - h(i))/dt;
    h(i) = integrate_S(h_start, S, i, dx);
    end
    %now recalcuate with pde
    S_guess = pde_Solve(S, Q, dt, dhdt);
end

%final loop 
for i = 1:nx
    dhdt(i) = integrate_S(h_start, S_guess, i, dx);
    h(i) = integrate_S(h_start, S_guess, i, dx); 
end
end 
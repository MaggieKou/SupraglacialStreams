function h = row_iterate(h_start, S_vec, dx, nx, dt, Q, err)
%start with previous row 
S = init_guess(S_vec, Q, dt);
h = zeros([1, length(S)]);
dhdt = zeros([1,length(S)]);
%start with integrating the slope along the line 
for i = 1:nx
    dhdt(i) = (integrate_S(h_start, S, i, dx) - h(i))/dt;
    h(i) = integrate_S(h_start, S, i, dx); % new h value 
end

%now recalcuate with pde
S_guess = pde_Solve(S_vec, Q, dt, dhdt); %solving from previous point

while (max(abs(S_guess - S)) > err)
    S = S_guess;
    for i = 1:nx
    dhdt(i) = (integrate_S(h_start, S, i, dx) - h(i))/dt;
    h(i) = integrate_S(h_start, S, i, dx);
    end
    %now recalcuate with pde
    S_guess = pde_Solve(S_vec, Q, dt, dhdt);
end

%final loop 
for i = 1:nx
    %dhdt(i) = integrate_S(h_start, S_guess, i, dx);
    h(i) = integrate_S(h_start, S_guess, i, dx); 
end
end 
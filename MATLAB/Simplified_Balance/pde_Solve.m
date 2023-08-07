function S_guess = pde_Solve(S_vec, Q, dt, dhdt)
S_guess = S_vec + dt.*(-5/2.*S_vec.^(12/5).*Q^(1/5) + 5/2.*Q^(-2/5).*S_vec.^(6/5).*dhdt); %guesses with assumption that h is the same a before
end
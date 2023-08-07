function S_guess  = init_guess(S_vec, h_diff, Q, dt)
%guesses current S based on previous value 
S_guess = S_vec + dt.*(-5/2.*S_vec.^(12/5).*Q^(1/5)); %guesses with assumption that h is the same a before
end
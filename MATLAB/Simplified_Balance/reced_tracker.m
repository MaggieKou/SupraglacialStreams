function [RHm, LHm, v_recR, v_recL, St_r, St_l, ht_r, ht_l] = reced_tracker(x, S, St, ht, slope_thresh, nt, nx, dt, dx)

RHm = zeros([1, nt]); %rhs migration rate
LHm = zeros([1, nt]); %lhs migration rate

St_l = zeros([1, nt]);
ht_l = zeros([1, nt]);
St_r = zeros([1, nt]);
ht_r = zeros([1, nt]);

for j = 1:nt
    for i = 1:nx
        if S(j,i) == slope_thresh %directly crosses 0 
            RHm(j) = x(i);
            St_r(j) = St(j,i);
            ht_r(j) = ht(j,i);
        else if i>1 && (S(j, i-1) - slope_thresh)*(S(j,i)-slope_thresh) <=0
            RHm(j) = x(i);  %just about crosses 0 at this point 
            St_r(j) = St(j,i);
            ht_r(j) = ht(j,i);
        end
        end
    end
end

for j = 1:nt
    for i = 1:nx
        if S(j,i) == slope_thresh %directly crosses
            LHm(j) = x(i);
            St_l(j) = St(j,i);
            ht_l(j) = ht(j,i);
            break
        else if i>1 && (S(j, i-1) - slope_thresh)*(S(j,i)-slope_thresh) <=0
            LHm(j) = x(i);
            St_l(j) = St(j,i);
            ht_l(j) = ht(j,i);
            break
        end
        end
    end
end

v_recR = zeros([1, nt]);
v_recL = zeros([1, nt]);

for j = 2:length(RHm)
    v_recR(j) = (RHm(j) - RHm(j-1))/dt;
    v_recL(j)= (LHm(j) - LHm(j-1))/dt;
end 
v_recR(1) = v_recR(2);
v_recL (1) = v_recL(2);


end

% Bingchen Liu Feb 25, 2025
% this code extract the vector with dx width from cxt_atx (21 by 21) matrix
% with +/- 1 dx 

function vec = extract_pm1(A,ind_row)
    vec = NaN(3,1);
    
    if ind_row-1>=1
        vec(1) = A(ind_row-1);
    end 

    vec(2) = A(ind_row);

    if ind_row+1 <= length(A)
        vec(3) = A(ind_row+1);
    end 
end 

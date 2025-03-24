% Bingchen Liu Feb 25, 2025
% this code extract the vector with dx width from cxt_atx (21 by 21) matrix
% with +/- 2 dx 

% function vec = extract_pm2(A,ind_row)
%     vec = NaN(5,1);
% 
%     for i = -2:2
%         if (ind_row+i>=1) && (ind_row+i<=length(A))
%             vec(i+3) = A(ind_row+i);
%         end 
%     end 
% end 


function vec = extract_pm2(A,ind_row)
    vec = NaN(9,1);

    for i = -4:4
        if (ind_row+i>=1) && (ind_row+i<=length(A))
            vec(i+5) = A(ind_row+i);
        end 
    end 
end 

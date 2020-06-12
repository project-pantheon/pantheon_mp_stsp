function [stops, c_req] = calculateStopsFromTreesIDs(trees_ID, n_r, n_c)
    %to expand

%     trees_IDs = [1 2 3 4 5 6 7 8 9]; 

    for i=trees_ID
        
%         stops(1) = mod(i,n_c) * (1+floor(i/n_c)) + n_r * floor(i/n_r);
%         stops(2) = mod(i,n_c) * (1+floor(i/n_c)) + n_r * floor(i/n_r) + 1;

        stops(1) = i + floor((i-0.01)/(n_c-1));%+  n_r * floor(i/n_r); %n° row
        stops(2) = i + floor((i-0.01)/(n_c-1)) +1;
        stops(3) = stops(1) + n_c; %ok
        stops(4) = stops(2) + n_c; %ok
    
    end

end


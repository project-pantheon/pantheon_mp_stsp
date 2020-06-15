function [required_vertex, c_req, associations] = calculateStopsFromTreesIDs(trees_IDs, n_c)
c_req = 0;
%     stops = zeros(4*length(trees_IDs));
    counter = 1;
    stops = [];
    
    associations = zeros(length(trees_IDs),5); %structured as [tree_id 4 stops]
     
    for i = trees_IDs

        first = i + floor((i-0.01)/(n_c-1));
        second =  i + floor((i-0.01)/(n_c-1)) +1;
        
        stops(end+1) = first;
        stops(end+1) = second;
        stops(end+1) = first + n_c; 
        stops(end+1) = second + n_c; 
        
        associations(counter,1) = i;
        associations(counter,2) = first;
        associations(counter,3) = second;
        associations(counter,4) = first + n_c;
        associations(counter,5) = second + n_c;
        
    
        counter = counter + 1;
    end
    
    required_vertex = unique(stops);
    
    c_req = histc(stops, required_vertex);
    
    

end


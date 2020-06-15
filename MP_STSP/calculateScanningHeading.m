function final_tour = calculateScanningHeading(tour, associations, stops_IDs, trees_IDs)

    %values are given relative to the heading of the robot
    % []
    
    
    final_tour = [tour, zeros(length(tour),4)];
    trees_to_scan = [];
    related_stops = [];
    
    for i=1:length(tour)
        
       if tour(i,2) == 1
           
          for j=1:size(associations,1)
          
              index = find(associations(j,2:end) == tour(i,1));

              if ~isempty(index)
                trees_to_scan = [trees_to_scan;associations(j,1)];
                related_stops = [related_stops;tour(i,1)];
              end
          
          end
         
       end
        
    end
    
    final_associations = [related_stops trees_to_scan ];
  
    
    stops_flags = zeros(length(tour),4);
    
    for i=1:size(final_associations,1)
        
        %simil atan2
        delta_x = trees_IDs(trees_to_scan(i),2) - stops_IDs(related_stops(i),2);
        delta_y = trees_IDs(trees_to_scan(i),3) - stops_IDs(related_stops(i),3);

        
        if delta_x > 0 && delta_y > 0 % [1 0 0 0]
              
            temp_mask =  [1 0 0 0];
            mask = rotate_mask(temp_mask,tour(tour(:,1) == related_stops(i) & tour(:,2) == 1,3));
            
        elseif delta_x < 0 && delta_y > 0 %[0 1 0 0]
            
             temp_mask =  [0 1 0 0];
          mask = rotate_mask(temp_mask,tour(tour(:,1) == related_stops(i) & tour(:,2) == 1,3));
        elseif delta_x < 0 && delta_y < 0 % [0 0 1 0]
            
             temp_mask =  [0 0 1 0];
              mask = rotate_mask(temp_mask,tour(tour(:,1) == related_stops(i) & tour(:,2) == 1,3));
            
        else
            
             temp_mask =  [0 0 0 1];
             
             mask = rotate_mask(temp_mask,tour(tour(:,1) == related_stops(i) & tour(:,2) == 1,3));
                     
        end
        
%        index = find(tour(:,1) == stops_IDs(i,1));
%        stops_flags(index,:) = stops_flags(index,:) + mask;
       
        % the second condition is necessary to distinguish between scanning
        % stop or simple passage to required stop
        final_tour(tour(:,1) == related_stops(i) & tour(:,2) == 1,4:7) = final_tour(tour(:,1) == related_stops(i)  & tour(:,2) == 1 ,4:7) + mask;
        
    end
    
    final_tour(:,3) = final_tour(:,3) *pi/2;
    
    %swap columns to have better layout and put -pi/2
    temp = final_tour(:,3);
    
    temp(temp > 4.7) = -pi/2;
    
    final_tour(:,3) = final_tour(:,2);
    final_tour(:,2) = temp;

end


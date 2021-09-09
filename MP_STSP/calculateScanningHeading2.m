function [final_tour, final_stops_IDs] = calculateScanningHeading2(tour, associations, stops_IDs, trees_IDs, delta_rows,delta_cols)

    %values are given relative to the heading of the robot
    % []
    
    
    %final_tour = [tour, zeros(size(tour,1),4)];
    trees_to_scan = [];
    related_stops = [];
    
    stops_IDs_heading=[];
    IDs_heading=size(stops_IDs,1);
    cnt=0;
    final_stops_IDs=stops_IDs;
    
    %check heading
    for i=1:size(tour,1)
        if(tour(i,2)==1 && ...                      %scan required
            (tour(i,3)==1 || tour(i,3)==3) )        %not feseable orientation
            if(i+1<=size(tour,1))                   %exists next waypoint(required)
                if(tour(i+1,3)==0 || tour(i+1,3)==2)%next feseable orientation
                    %build 2 fake waypoints in the opposite direction                    
                    if tour(i,3)==1 && tour(i+1,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==1 && tour(i+1,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i+1,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i+1,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    end
                    %use orientation of next waypoint
                    tour(i,3)=tour(i+1,3);
                    
                elseif(i+2<=size(tour,1) && (tour(i+2,3)==0 || tour(i+2,3)==2)) %next next position
                    %build 2 fake waypoints in the same direction
                     if tour(i,3)==1 && tour(i+2,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==1 && tour(i+2,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i+2,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i+2,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    end
                    
                    %use orientation of next waypoint
                    tour(i,3)=tour(i+2,3);
                    
                else
                    %build 2 fake waypoints in the same direction
                    if tour(i,3)==1 && tour(i-1,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==1 && tour(i-1,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i-1,3)==0
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    elseif tour(i,3)==3 && tour(i-1,3)==2
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                        cnt=cnt+1;
                        stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                    end
                    
                    %use opposite direction of prev waypoint
                    tour(i,3)=4-tour(i-1,3);
                    
                end
            else
                %build 2 fake waypoints in the same direction
                if tour(i,3)==1 && tour(i-1,3)==0
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                elseif tour(i,3)==1 && tour(i-1,3)==2
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows*2/3 tour(i,1) tour(i,3)];
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)-delta_rows/3 tour(i,1) tour(i,3)];
                elseif tour(i,3)==3 && tour(i-1,3)==0
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)+delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                elseif tour(i,3)==3 && tour(i-1,3)==2
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows*2/3 tour(i,1) tour(i,3)];
                    cnt=cnt+1;
                    stops_IDs_heading=[stops_IDs_heading; IDs_heading+cnt stops_IDs(tour(i,1),2)-delta_cols stops_IDs(tour(i,1),3)+delta_rows/3 tour(i,1) tour(i,3)];
                end
                
                %use opposite direction of prev waypoint
                tour(i,3)=4-tour(i-1,3);
            end
        end
    end
    if ~isempty(stops_IDs_heading)
        final_stops_IDs=[final_stops_IDs; stops_IDs_heading(:,1:3)];
    end
    final_tour = [tour, zeros(size(tour,1),4)];
    
    for i=1:size(tour,1)
        
       if tour(i,2) == 1
           
          for j=1:size(associations,1)
          
              index = find(associations(j,2:end) == tour(i,1));

              if ~isempty(index)
                trees_to_scan = [trees_to_scan;associations(j,1)-1];
                related_stops = [related_stops;tour(i,1)];
              end
          
          end
         
       end
        
    end
    
    final_associations = [related_stops trees_to_scan ];
    
    stops_flags = zeros(size(tour,1),4);
    
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
    
    %add extra waypoint
    final_tour_real = zeros(size(tour,1),7);
    j=1;
    for i=1:size(final_tour,1)
        if j<size(stops_IDs_heading,1)
            if final_tour(i,1)== stops_IDs_heading(j+1,4)
                final_tour_real(i+j-1,:)=[stops_IDs_heading(j,1) 0 stops_IDs_heading(j,5) 0 0 0 0];
                j=j+1;
                final_tour_real(i+j-1,:)=[stops_IDs_heading(j,1) 0 stops_IDs_heading(j,5) 0 0 0 0];
                j=j+1;
            end
        end
        final_tour_real(i+j-1,:)=final_tour(i,:);
    end
    final_tour=final_tour_real;
    
    final_tour(:,3) = final_tour(:,3) *pi/2;
    
    %swap columns to have better layout and put -pi/2
    temp = final_tour(:,3);
    
    temp(temp > 4.7) = -pi/2;
    
    final_tour(:,3) = final_tour(:,2);
    final_tour(:,2) = temp;

end

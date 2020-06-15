function [tour_robot_heading] = findRobotHeading(tour, n_cols)

    %heading robot can have 4 values:
    % --> 0
    % ^ 1
    % <-- 2
    % v 3
    tour_robot_heading = [tour, zeros(length(tour),1)];
    
    diff_ids = diff(tour(:,1));
    
    for i=2:size(tour,1)
       
        if diff_ids(i-1) == 1
            
            tour_robot_heading(i,3) = 0;
            
        elseif diff_ids(i-1) == -1
            
            tour_robot_heading(i,3) = 2;
            
        elseif diff_ids(i-1) == n_cols
            
            tour_robot_heading(i,3) = 1;
        
        else
            
            tour_robot_heading(i,3) = 3;
            
        end
        
    end
    
    tour_robot_heading(1,3) = tour_robot_heading(2,3);



end


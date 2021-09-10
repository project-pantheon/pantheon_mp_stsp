function [edges_list, dist] = calculateEdgesList(stops_coord,nStops,grid)

    threshold=0.3;% m

    X = stops_coord(:,1);
    Y = stops_coord(:,2);
    
    edges_list = nchoosek(1:nStops,2); %COMPLETE GRAPH


    edges_list = [edges_list; edges_list(:,2), edges_list(:,1)];
    
    dist = hypot(Y(edges_list(:,1)) - Y(edges_list(:,2)), ...
                 X(edges_list(:,1)) - X(edges_list(:,2)));
             
            

    id1 = dist > max(grid(1)+threshold,grid(2)+threshold);
    dist(id1) = [];
    edges_list(id1,:) = [];
   


end


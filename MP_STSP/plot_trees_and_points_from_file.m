function [XY,T, origin_utm] = plot_trees_and_points_from_file(map_file,orientation,nRobots,n_r,n_c,d_r,d_c,required_vertex)

    close all;
    
    if n_c < 10 && n_r < 10
        trees_size = 100;
        stops_size = 100;
    else
        trees_size = 20;
        stops_size = 20;
        
    end
    
    
    val_tab_map = readtable(map_file);
    val_map=zeros(size(val_tab_map,1),2);
    val_map_utm=zeros(size(val_tab_map,1),2);
    val_map(:,1)=val_tab_map.Var2;
    val_map(:,2)=val_tab_map.Var3;
    origin=val_map(1,:);
    origin_utm=ll2utm(origin);
    for i=1:size(val_tab_map,1)
        temp=(rotz(rad2deg(-orientation))*[ll2utm(val_map(i,:))-origin_utm 0]')';
        val_map_utm(i,:)=temp(1:2)';
    end
    
    row_mask = 0:d_c:(n_c-1)*d_c;
    col_mask = 0:d_r:(n_r-1)*d_r;
    [x_stops, y_stops] = meshgrid(row_mask,col_mask);

    XY = [reshape(x_stops',[n_r*n_c 1]) reshape(y_stops',[n_r*n_c 1])];

    row_mask = 0:d_c:(n_c-2)*d_c;
    col_mask = 0:d_r:(n_r-2)*d_r;
    [x_trees, y_trees] = meshgrid(row_mask,col_mask);

    x_trees = x_trees + d_c/2; 
    y_trees = y_trees + d_r/2; 
    T = [reshape(x_trees',[(n_r-1)*(n_c-1) 1]) reshape(y_trees',[(n_r-1)*(n_c-1) 1])];
    
    T_utm=zeros(n_c-1,n_r-1);
    for i=1:n_c-1
        for j=1:n_r-1
            x=[val_map_utm(i+n_c*(j-1),1) val_map_utm(i+n_c*(j-1)+1,1) val_map_utm(i+n_c*j,1) val_map_utm(i+n_c*j+1,1)];
            y=[val_map_utm(i+n_c*(j-1),2) val_map_utm(i+n_c*(j-1)+1,2) val_map_utm(i+n_c*j,2) val_map_utm(i+n_c*j+1,2)];
            T_utm(i,1) = mean(x);
            T_utm(i,2) = mean(y);
        end
    end

    %virtual values overwritten
    XY=val_map_utm;
    T=T_utm;
    
    depot = [-5 -5];
    XY = [depot;XY];
     
    for kk=1:nRobots+1 %why is n+1 here?
        
        figure(kk)  
        hold on
        title(['robot n°', num2str(kk)])
        % Add the trees to the map
        scatter(T(:,1),T(:,2),trees_size,'o','filled','g');
        % Add the stops to the map
        scatter(XY(:,1),XY(:,2),stops_size,'o','filled','b')

        for j=required_vertex
            scatter(XY(j,1),XY(j,2),stops_size*1.2,'o','filled','r')
        end
          
        if n_r < 10 && n_c < 10
            for i=1:length(XY)
                text(XY(i,1) + 0.3,XY(i,2) + 0.35,int2str(i),'Color', 'b')
            end

            for i=1:length(T)

                text(T(i,1) + 0.3,T(i,2) + 0.35,int2str(i),'Color', 'g')

            end
        end
    
    end


end


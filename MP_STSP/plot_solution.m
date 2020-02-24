function  [fig,sel_edges] = plot_solution(fig,sol,idxs,stopsLat,stopsLon,nRobots, ...
         n_rows,n_cols,delta_rows,delta_cols,required_vertex)


    
    single_tour_len = length(sol)/nRobots;

    for i=1:nRobots
        
        x_mtsp(i,:) = sol((i-1)*single_tour_len +1:i*single_tour_len) ;
        sel_edges{i} = round(idxs(sol >= .99,:));
        segments{i} = find(round(x_mtsp(i,:))); 
    
        Lat{i} = zeros(3*length(segments{i}),1);
        Lon{i} = zeros(3*length(segments{i}),1);
    
        for ii = 1:length(segments{i})
            start = idxs(segments{i}(ii),1);
            stop = idxs(segments{i}(ii),2);

            % Separate data points with NaN's to plot separate line segments
            Lat{i}(3*ii-2:3*ii) = [stopsLat(start); stopsLat(stop); NaN];
            Lon{i}(3*ii-2:3*ii) = [stopsLon(start); stopsLon(stop); NaN];  
        end
    end
    
    colors = 'cmyk';

    for i=1:nRobots
        figure;
       plot_trees_and_points(n_rows,n_cols,delta_rows,delta_cols,required_vertex);

        plot(Lat{i},Lon{i},'Color',colors(i),'LineWidth',2);

    end

end


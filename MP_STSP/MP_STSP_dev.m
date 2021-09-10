clear all;
close all;
clc;

%format long
addpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux')

solution_found = 0;
solution_errors = 0;
mipgap_=.001;
offset_=0;

while ~solution_found && ~solution_errors
    
    try
        
        %% Parse input file 
        
        %Map
        fnamemap = 'input/maps_config.json';
        valmap = jsondecode(fileread(fnamemap));
        map_index=0;
        
        %Mission
        fname = 'input/mission.json';
        val = jsondecode(fileread(fname));
        
        i=1;
        while map_index==0 && i<=size(valmap,1)
            if strcmp(valmap(i).id,val.maps_config_id)
                map_index=i;
                id_mission_base=val.id_mission_base;
            end
            i=i+1;
        end
        
        %Map
        map_file = "input/"+ valmap(map_index).map_file;
        orientation = valmap(map_index).field.orientation;
        
        %Trees
        fname = "input/"+ valmap(map_index).mapping;
        val_tab_trees = readtable(fname);
        val_trees.id=val_tab_trees.Var1;
        val_trees.name=[""];
        i=1;
        for i=1:size(val_tab_trees.Var2(:),1)
            val_trees.name(i)=val_tab_trees.Var2{i};
        end
        
        n_rows_trees = valmap(map_index).field.trees_rows;
        n_cols_trees = valmap(map_index).field.trees_cols;
        delta_rows = valmap(map_index).field.delta_rows;
        delta_cols = valmap(map_index).field.delta_cols;
        
        delta_cols = delta_cols - offset_;
        
        full_trees=val.full_trees; % if full mode is set, all ids will be ignored
        
        if full_trees
            trees_IDS_to_scan = linspace(1, n_rows_trees * n_cols_trees, n_rows_trees * n_cols_trees );
        else
            trees_IDS_to_scan = val_trees.id(ismember(val_trees.name,val.targets))';
        end
        nRobots = val.nrobots;
        
        
        %%
        
        which_solver = 2; %1 for intlinprog, 2 for cplexmilp
        
        %plot parameters
        
        colors = 'kkkkk';
        
        %grid parameters
        n_rows = n_rows_trees + 1;
        n_cols = n_cols_trees + 1;
        
        
        grid = [delta_rows, delta_cols];
        
        %trees_IDS_to_scan = [1 6 5 10];
        
        [required_vertex, c_req, associations] = calculateStopsFromTreesIDs(trees_IDS_to_scan, n_cols);
        
        Tmax_cost = 10; %if we push too much the time, the robots tend to split exactly the number of required nodes
        
        if nRobots == 1 %single robot Tcapture has no role
            
            Tmax_cost = 1;
            
        end
        
        c1 = 40; %single stop cost
        c2 = 5; %traverse edge cost
        
        c_req = 40*[c_req];
        
        %fix ids
        required_vertex = 1+[required_vertex];
        associations=1+associations;
        
        
        nStops = n_rows * n_cols;
        X = zeros(nStops,1);
        Y = X;
        
        %% Plot
        
        required_vertex = unique(required_vertex); %eliminates repetitions (if any)
        %depot_indices = find(required_vertex == 1);
        %required_vertex(depot_indices) = [];%take out depot if present
        %c_req(depot_indices) = [];
        
        nRequired = length(required_vertex);
                
        %[XY,T] = plot_trees_and_points(nRobots,n_rows,n_cols,delta_rows,delta_cols,required_vertex);
        [XY,T, origin_utm] = plot_trees_and_points_from_file(map_file,orientation,nRobots,n_rows,n_cols,delta_rows,delta_cols,required_vertex);
        
        %%
        
        %get single parameters vectors
        [edges_list, dist] = calculateEdgesList(XY,nStops+1,grid);
        %edges_list_r
        edges_len = length(edges_list);
        
        %add depot
        depot_edges=[1 2];
        rev_dep_edges=[2 1];
        depot_cost = zeros(size(depot_edges,1),1);
        dist = [depot_cost ;dist(1:edges_len/2,:);depot_cost; dist(edges_len/2+1:end,:)];
        edges_list = [depot_edges ;edges_list(1:edges_len/2,:); rev_dep_edges ;edges_list(edges_len/2+1:end,:)];
        
        
        
        [delta_plus , delta_minus, Fa] = calculateFlowVariables(nStops+1,edges_list);
        
        V = calculateSplitVariables(nRequired);
        
        XY_with_IDs = [linspace(1,length(XY),length(XY))' XY];
        T_with_IDs = [linspace(1,length(T),length(T))' T];
        
        
        %increase variable for nRobots
        
        if nRobots > nRequired
            
            fprintf('Too many robots for this application.\n');
            fprintf('Reducing the number of robots.');
            
            while nRobots > nRequired
                
                nRobots = nRobots - 1;
                
            end
            
        end
        
        edges_list = repmat(edges_list,nRobots,1);
        
        edges_len = length(edges_list);
        
        dist = repmat(dist,nRobots,1);
        
        Fa = repmat(Fa,nRobots,1);
        
        flow_len = length(Fa);
        
        V = repmat(V,nRobots,1);
        
        v_len = length(V);
        
        solution_len = edges_len + flow_len + v_len + 1; % last is Tmax
        single_rob_len = edges_len/nRobots;
        
        %ILP variables and cost function
        intcon = 1:solution_len;
        
        %CONSIDER ALSO PATH MINIMIZATION (not working otherwise)
        % cost_fun = [ dist ; Fa ; V ; Tmax_cost];
        
        % cost_fun = [ dist ; zeros(length(Fa)+length(V),1) ; Tmax_cost];
        
        
        % SPARSE VECTOR
        cost_fun = spalloc(length(dist)+length(Fa)+length(V)+1,1,1);
        %cost_fun = [ dist; zeros(length(Fa)+length(V),1); Tmax_cost];
        cost_fun(end) = Tmax_cost;
        
        
        %% mSTSP additional variables
        
        [Aeq_v , beq_v, A_v, b_v] = calculateVCnstrs(solution_len,single_rob_len,v_len,nRobots,nRequired);
        %OK PARAMETRIZED
        
        [A_t, b_t] = calculateTimeCnstrs(nRobots,c_req,c2,solution_len,single_rob_len,nRequired);
        %OK PARAMETRIZED
        
        %% Equality Constraints
        
        [Aeq , beq] = calculateEqCnstrs(nRobots,nStops,edges_len,delta_plus,delta_minus,nRequired,required_vertex);
        %OK PARAMETRIZED
        
        %% Not equality constraints
        
        [A , b] = calculateNEqCnstrs(nRobots,nStops,edges_len,delta_plus,required_vertex);
        
        
        %% Lower Bounds (upper bounds already considered)
        
        lb = spalloc(length(intcon),1,0);
        
        % lb = zeros(length(intcon),1);
        
        
        %% Solve Optimization
        
        %impose flow at depot = nRequired (SHOULDNT BE REQUIRED-->DEBUG OTHERS)
        % Aeq_d = spalloc(1,solution_len,nRobots*4) ;
        %
        % Aeq_d(1,edges_len+1:edges_len+edges_len/2) = delta_plus{1};
        % Aeq_d(1,edges_len+edges_len/2+1:edges_len+edges_len) = delta_plus{1};
        %
        %
        % beq_d = nRequired;
        
        
        Aeq = [Aeq;Aeq_v];%;Aeq_d];
        beq = [beq;beq_v];%;beq_d];
        
        A = [A;A_t];%A_v;
        b = [b;b_t];%b_v;
        
        if which_solver == 1
            
            opts = optimoptions('intlinprog','Display','iter');
            [x_tsp,costopt,exitflag,output] = intlinprog(cost_fun,intcon,A,b,Aeq,beq,lb,[],opts);
            
        else
            %CPLEX
            options = cplexoptimset('cplex');
            options.display = 'on';
            options.mip.limits.nodes = 100000;
            %     options.mip.strategy.search = 3;
            %     options.mip.strategy.search = 1; %1 branch&cut ---> search other approaches
            options.mip.tolerances.absmipgap = 1;
            options.mip.tolerances.mipgap = mipgap_;
            options.mip
            
            
            
            ctype = 0;
            
            for i=1:edges_len
                
                ctype = strcat(ctype,'B');
                
            end
            for i=edges_len+1:2*edges_len
                
                ctype = strcat(ctype,'I');
                
            end
            for i=2*edges_len+1:2*edges_len+v_len
                
                ctype = strcat(ctype,'B');
                
            end
            
            ctype = strcat(ctype,'I');
            
            %OPTIMIZATION POSSIBLES: SOSvariables (i vri che sono mutualmente
            %esclusivi)
            %UPPER BOUNDS PER TUTTI
            pause(5)
            [x_tsp, fval, exitflag, output] = cplexmilp(cost_fun, A, b, Aeq, beq, [],[],[], lb, [], ctype,[], options);
            
        end
        
        
        %% write solution
        % vir = round(x_tsp(nRobots*single_rob_len+1:end-1));
        
        segments = find(x_tsp); % Get indices of lines on optimal path
        lh = zeros(nStops,1); % Use to store handles to lines on plot
        
        sol_edges = cell(nRobots,1);
        sol_flows = cell(nRobots,1);
        missions_time = zeros(nRobots);
        
        for i=1:nRobots
            
            sol_edges{i} = round(find_edges(x_tsp((i-1)*single_rob_len +1:i*single_rob_len),edges_list));
            sol_flows{i} = round(find_flows(x_tsp(edges_len +(i-1)*single_rob_len +1:edges_len +i*single_rob_len),edges_list));
            sol_v(i,:) = round(x_tsp(2*edges_len+(i-1)*v_len/nRobots+1:2*edges_len+i*v_len/nRobots));
            
            
            Vri = required_vertex(sol_v(i,:) == 1);
            updatePlot(i,x_tsp((i-1)*edges_len/nRobots +1:i*edges_len/nRobots),edges_list(1:edges_len/nRobots,:),XY,colors(i),Vri);
            
            
            tour{i} = findTour(sol_edges{i},required_vertex,Vri);
            tour_real{i} = findTourReal(sol_edges{i},required_vertex,Vri);
            
            %T_tot(i) = calculateFinalTime(x_tsp((i-1)*edges_len/nRobots +1:i*edges_len/nRobots),Vri,required_vertex,c_req);
            T_tot(i) = calculateFinalTime(tour{i},Vri,required_vertex,c_req);
            fprintf(2,'Total time for robot # %i\n',i);
            fprintf(2,'%i [sec]\n' , T_tot(i));
            
        end
        
        
        % Robots heading are given as they arrive to the stop
        % Heading at the depot (stop n°1) has been left as unknown. It has to be given
        % as input param. For now it is equal to the first stop the robots do.
        XY_with_IDs_temp=XY_with_IDs;
        for i=1:nRobots
            
            tour_with_robot_heading{i} = findRobotHeading(tour{i},n_cols);
            %final_tour{i} = calculateScanningHeading(tour_with_robot_heading{i},associations,XY_with_IDs, T_with_IDs);
            [final_tour{i}, final_stops_IDs{i}] = calculateScanningHeading2(tour_with_robot_heading{i},associations,XY_with_IDs_temp, T_with_IDs,delta_rows,delta_cols);
            new_waypoints=size(final_stops_IDs{i},1)-size(XY_with_IDs_temp,1);
            XY_with_IDs_temp=final_stops_IDs{i};
            figure(i)
            scatter(final_stops_IDs{i}(size(XY_with_IDs_temp,1)-new_waypoints+1:end,2),final_stops_IDs{i}(size(XY_with_IDs_temp,1)-new_waypoints+1:end,3),100,'o','filled','g')
        end
        
        
        for i=1:nRobots
            for j=1:size(final_stops_IDs{i},1)
                temp=rotz(rad2deg(orientation))*[final_stops_IDs{i}(j,2:3) 0]'+[origin_utm'; 0];
                final_stops_IDs{i}(j,2:3)= utm2ll(temp(1),temp(2),33);
            end
        end
        
        
        
        %% Write solution files
        
        fileID = fopen('edges.txt', 'w');
        for i=1:length(edges_list)
            fprintf(fileID,'(%d, %d),\n', edges_list(i,1)-1,edges_list(i,2)-1);
        end
        fclose(fileID);
        
        % Final tour is composed as:
        % 1 column: path consisting into all the stops the i-th robot has to do
        % 2 column: heading of the robot for each single stop (discretized into 4)
        % 3 column: boolean vector (0 scanning not required, 1 scanning required)
        % 4:7 columns: scanning positions relative to robot heading of column 2
        
        for i=1:nRobots
            writematrix(cell2mat(final_tour(i)), "/home/data/mongodb/robots/UGV/missions/"+id_mission_base+"_"+num2str(i)+"/tour.csv");
            writematrix(cell2mat(final_stops_IDs(i)), "/home/data/mongodb/robots/UGV/missions/"+id_mission_base+"_"+num2str(i)+"/map_processed.csv");
        end
        
        solution_found = 1
    catch
        %if solver generates some errors, the mingap will be scale by 10%
        %or it will be introduced an offset on one delta
        if offset_~= 0
            offset_= 0
            mipgap_=mipgap_*0.9

            if (mipgap_ < 0.01)
                solution_errors = 1
            end
        else
            offset_=.0001
        end
    end
    
end


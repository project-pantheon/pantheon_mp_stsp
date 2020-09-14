clear all;
close all;
clc;

%% Parse input file

%fname = 'Mission.json';
fname = 'full_orchard_30x20.json';
val = jsondecode(fileread(fname));


n_rows_trees = str2double( table2array( struct2table(val.field{1}))); 
n_cols_trees = str2double( table2array(struct2table(val.field{2}))); 
delta_rows = str2double( table2array(struct2table(val.field{3}))); 
delta_cols = str2double( table2array(struct2table(val.field{4})));

if contains(fname, 'full_orchard')
    trees_IDS_to_scan = linspace(1, n_rows_trees * n_cols_trees, n_rows_trees * n_cols_trees );
else
    trees_IDS_to_scan = str2double(struct2cell(val.trees));
end
nRobots = str2double(val.nrobots);



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

nStops = n_rows * n_cols;
X = zeros(nStops,1); 
Y = X;

%% Plot

required_vertex = unique(required_vertex); %eliminates repetitions (if any)
depot_indices = find(required_vertex == 1);
required_vertex(depot_indices) = [];%take out depot if present
c_req(depot_indices) = [];

nRequired = length(required_vertex);

[XY,T] = plot_trees_and_points(nRobots,n_rows,n_cols,delta_rows,delta_cols,required_vertex);

%% 

%get single parameters vectors
[edges_list, dist] = calculateEdgesList(XY,nStops,grid);

[delta_plus , delta_minus, Fa] = calculateFlowVariables(nStops,edges_list);

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
%     options.mip.tolerances.absmipgap = 1;
%     options.mip.tolerances.mipgap = 0.001;
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

    T_tot(i) = calculateFinalTime(x_tsp((i-1)*edges_len/nRobots +1:i*edges_len/nRobots),Vri,required_vertex,c_req);
    fprintf(2,'Total time for robot # %i\n',i); 
    fprintf(2,'%i [sec]\n' , T_tot(i));
    
end


% Robots heading are given as they arrive to the stop
% Heading at the depot (stop n°1) has been left as unknown. It has to be given
% as input param. For now it is equal to the first stop the robots do.

for i=1:nRobots

    tour_with_robot_heading{i} = findRobotHeading(tour{i},n_cols);
    final_tour{i} = calculateScanningHeading(tour_with_robot_heading{i},associations,XY_with_IDs, T_with_IDs);

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

 writematrix(cell2mat(final_tour(i)), num2str(i) + "_robot_tour.csv");

end




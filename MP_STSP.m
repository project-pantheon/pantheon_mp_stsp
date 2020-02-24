clear all;
close all;
%clc;


%% USED PARAMETERS (17-02)
% orchard 9x7
% planting patter 5x5
% required_vertex = [10 11 17 18 12 19 26 27 33 34 44 45 51 52 48 49 55 56 62 63]; 
% nRobots from 1 to 3 (4 is too much..)
% cplex tolerance 0.5% convergence
% single stop c =
% c'è ridondanza....

which_solver = 2; %1 for intlinprog, 2 for cplexmilp 

%plot parameters

colors = 'kkkkk';

%grid parameters

n_rows = 7; % 3 
n_cols = 9; % 4

delta_rows = 5;
delta_cols = 5;

grid = [delta_rows, delta_cols];

nRobots = 1; %n° robots

Tmax_cost = 10; %if we push too much the time, the robots tend to split exactly the number of required nodes

if nRobots == 1 %single robot Tcapture has no role
    
    Tmax_cost = 1;
    
end

c1 = 40; %single stop cost
c2 = 5; %traverse edge cost

nStops = n_rows * n_cols;
X = zeros(nStops,1); 
Y = X;

% to have all stops (TSP-like)
% required_vertex = 2:nStops;

if n_rows == 2
    
    required_vertex = [2 3];
    
elseif n_rows == 3
    
    required_vertex = [6 7 8];
    
elseif n_rows == 4
    
%   required_vertex = [2 5 7 8 9 10 12];
%   required_vertex = [8 9 10 12 13 14]; %OK!!
    required_vertex = [2 6 8 9 10 4 13 14 16];
    
elseif n_rows == 5
    
    required_vertex = [3 4 8 9 14 15 19 20 25 24];

    c_req = [1 1 1 1 1 1 2 2 1 1];

  
elseif n_rows == 6
    
    required_vertex = [5 6 11 12 26 27 20 21 30 29 35 36];
    
elseif n_rows == 9 || n_cols == 9
    
%     required_vertex = [31 32 35 36 40 41 44 45 62 63 64 65 71 72 73 74 80 81];
%       c_req =           [1  1  1  1  1  1  1  1  1  1  1  1  2  2  1  1  1  1];
    required_vertex = [22 23 30 31 32 39 40 43 44 45 52 53 54 61 62 63]; 
    c_req =     c1 *  [1  1  1  2  1  1  1  1  2  1  2  4  2  1  2  1];
    
elseif n_rows == 10
    
%     31--> 1
%     32--> 1
%     41--> 1
%     42--> 1
%     35--> 1
%     36--> 1
%     45--> 1
%     46--> 1
%     47--> 1
%     48--> 1
%     57--> 1
%     58--> 1
%     62--> 1
%     63--> 2
%     64--> 1
%     72--> 1
%     73--> 2
%     74--> 1
%     79--> 1
%     89--> 2
%     80--> 1
%     90--> 2
%     99--> 1
%     100--> 1
    
    
    required_vertex = [31 32 41 42 35 36 45 46 62 63 64 72 73 74 79 89 80 90 99 100];

    c_req =           [1  1  1  1  1  1  1  1  1  2  1  1  2  1  1  2  1  2  1  1];
    
else
    
     required_vertex = [ 24 3 25 7 14 40 13 12 5 28 6 15 49 58 80 75 99 102];
 
end

required_vertex = unique(required_vertex); %eliminates repetitions (if any)
nRequired = length(required_vertex);

[XY,T] = plot_trees_and_points(nRobots,n_rows,n_cols,delta_rows,delta_cols,required_vertex);

%get single parameters vectors
[edges_list, dist] = calculateEdgesList(XY,nStops,grid);

[delta_plus , delta_minus, Fa] = calculateFlowVariables(nStops,edges_list);

V = calculateSplitVariables(nRequired);


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


%% 

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

%     missions_time(i) = c1*sum(sol_v(i,:)) + c2*length(sol_edges{i});
%     fprintf('Total time for robot # %i\n',i); 
%     fprintf('%i [sec]\n' , missions_time(i));
%     
    
    tour{i} = findTour(sol_edges{i},required_vertex,Vri);

    T_tot(i) = calculateFinalTime(x_tsp((i-1)*edges_len/nRobots +1:i*edges_len/nRobots),Vri,required_vertex,c_req);
    fprintf(2,'Total time for robot # %i\n',i); 
    fprintf(2,'%i [sec]\n' , T_tot(i));
    
end

%fprintf( 2,'T\n' )
%fprintf(2,'%i [sec]\n', x_tsp(end));



fileID = fopen('edges.txt', 'w');
for i=1:length(edges_list)
    fprintf(fileID,'(%d, %d),\n', edges_list(i,1)-1,edges_list(i,2)-1);
end
fclose(fileID);




function [A , b] = calculateNEqCnstrs(nRobots,nStops,edges_len,delta_plus,required_vertex)

    sol_len = 2*edges_len+nRobots*length(required_vertex) +1; %it considers already the nRobots
    robot_len = edges_len/nRobots;
    nRequired = length(required_vertex)+1; %we are including depot now

    % -xe <= -1 for i in Vr including depot
    
    ones_in_A_req = nRobots*( sum(sum([delta_plus{[1 required_vertex]}])));

    A_req = spalloc(nRequired,sol_len,ones_in_A_req) ;

    k = 1;
    
    for i=[1 required_vertex]
            
        for j=1:nRobots

            A_req(k,(j-1)*robot_len+1:j*robot_len) = -delta_plus{i};
            
%             A_req(k,edges_len/2 + 1:edges_len/2 + robot_len) = -delta_plus{i};
        end
        
        k = k+1;
    
    end
    

    b_req = -ones(nRequired,1);
    b_req(1) = -nRobots; %depot must be visited more than nRobots times
    
    % fe - (nr-1)xe <= 0

    A_f_bound = spalloc(edges_len,sol_len,2*edges_len) ;


    for i=1:edges_len
    
        A_f_bound(i,i) = -(nRequired); %depots
        A_f_bound(i,i + edges_len) = 1;
            
    end

   
    b_f_bound = spalloc(edges_len,1,0);
%     b_f_bound = zeros(edges_len,1);


    A = [A_req;A_f_bound];

    A_to_show = full(A);

    b = [b_req;b_f_bound];




end


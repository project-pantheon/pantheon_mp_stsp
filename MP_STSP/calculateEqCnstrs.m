function [Aeq , beq] = calculateEqCnstrs(nRobots,nStops,edges_len,delta_plus,delta_minus,nRequired,required_vertex)

    %structure of A: |A_robot1 edges |
    %                |A_robot2 edges |
    %                       ..
    %                |  A_robot1 fa  |
    %                |  A_robot2 fa  |
    %                       ..
    
    %all lengths
    sol_len = 2*edges_len+nRobots*nRequired + 1; %it considers already the nRobots
    robot_len = edges_len/nRobots;
    v_len = nRobots*nRequired;
   
        
    % sum(xe) = sum(xe) for each robot %OK
    
    %TO PARAMETRIZE FOR nRobots

    Aeq_x = spalloc(nRobots*nStops,sol_len,2*edges_len);
    
    for j=1:nRobots
        
        for i=1:nStops

            Aeq_x(i+(j-1)*nStops,(j-1)*robot_len+1:j*robot_len) = ( delta_plus{i} - delta_minus{i} )';

        end
        
    end
    
    beq_x = spalloc(nRobots*nStops,1,0);

%     beq_x = zeros(nRobots*nStops,1);

    %flows constraints
    %fe for i in Vr\{1} and fe in V\Vr
    
    %TO PARAMETRIZE

    
    Aeq_f = spalloc(nRobots*(nStops-1),sol_len,2*edges_len + nRobots*nRequired-nRobots*4); %depot will always have 2+2 edges (thus the 4)
    
    for j=1:nRobots
        
        for i=1:nStops-1
        
        Aeq_f(i+(j-1)*(nStops-1),edges_len+(j-1)*robot_len+1:edges_len+j*robot_len) = ( delta_minus{i+1} - delta_plus{i+1} )'; %avoids depot (at 1)
        
        end
        
    end
    
     
    kk = 0; %counter that associates the i^th required with the i^th v
    for j=1:nRobots
        
        for i=required_vertex
            
            Aeq_f(i-1 +(j-1)*(nStops-1), 2*nRobots*robot_len+1+kk) = -1;
            kk = kk+1;
        
        end
        
    end
    
    beq_f = spalloc(nRobots*(nStops-1) ,1 , 0);

%     beq_f = zeros( nRobots*(nStops-1) ,1);
    
    Aeq = [Aeq_x;Aeq_f];

    beq = [beq_x; beq_f];
    
        
    %OK UNTIL HERE


end


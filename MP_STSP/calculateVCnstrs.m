function [Aeq_v , beq_v, A_v , b_v] = calculateVCnstrs(sol_len,single_rob_len,v_len,nRobots,nRequired)

    %equality constraints
    %sum of all the vir = |Vr|
    
    %WATCH OUT, HERE I'M NOT CONSIDERING DEPOT 1 FOR EVERY ROBOT AS REQUIRED
    
    Aeq_v = spalloc(1,sol_len,nRobots*nRequired);
    
    Aeq_v(1,2*single_rob_len*nRobots+1:2*single_rob_len*nRobots+nRobots*nRequired) = ones(1,nRobots*nRequired);
    
    beq_v = nRequired;%(nRobots -1) is necessary to consider the depot for every robot
    
  
    %upper bound for vir to have binary values (TO INCLUDE IN CPLEX ctype
    %VECTOR)
    
    A_v = spalloc(nRobots*nRequired,sol_len,nRobots*nRequired);
    
    for i=1:nRobots*nRequired
       
        A_v(i,2*single_rob_len*nRobots+i) = 1;
        
    end
    
    
    b_v = ones(nRobots*nRequired,1);
    
    % sum 2-by2 vir = 1 (mutually exclusive) except (EVENTUALLY NEXT FORMULATION) sum depot =nRobots
    %ALSO THIS CAN BE INCLUDED IN CPLEX FORMULATION (WITH SOSvariables)
    
    Aeq_e = spalloc(nRequired,sol_len,nRobots*nRequired);
    
    for i=1:nRequired
    
        for j=1:nRobots 

            Aeq_e(i,2*single_rob_len*nRobots+i+(j-1)*v_len/nRobots) = 1;


        end
    
    end
    
    
    beq_e = ones(nRequired,1);
    
     %Aeq_v = [Aeq_v ; Aeq_e];
Aeq_v = Aeq_e;
     %beq_v = [beq_v ; beq_e];
beq_v = beq_e;
    
end


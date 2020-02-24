function [A_t, b_t] = calculateTimeCnstrs(nRobots,c1, c2,sol_len,single_rob_len,nRequired)

    %for each robot we define a tmax
    A_t = spalloc(nRobots ,sol_len,(single_rob_len+nRequired +1) * nRobots); %added ub on Tmax

    for k=1:nRobots
        
        if length(c1) > 1
        
            A_t(k,2*single_rob_len*nRobots+1+(k-1)*nRequired:2*single_rob_len*nRobots+k*nRequired) = c1.*ones(1,nRequired);
       
        else
            
            A_t(k,2*single_rob_len*nRobots+1+(k-1)*nRequired:2*single_rob_len*nRobots+k*nRequired) = c1*ones(1,nRequired);
            
        end
        
       %let's avoid the edges cost repetition
       A_t(k,(k-1)*single_rob_len+1:k*single_rob_len) = c2*ones(1,single_rob_len);
        
       A_t(k,end) = -1;
        
    end

    b_t = zeros(nRobots,1); 
    
    
    %with upper bound on Tmax:
%     A_t = spalloc(nRobots +1 ,sol_len,(single_rob_len+nRequired +1) * nRobots); %added ub on Tmax
% 
%     for k=1:nRobots
%         
%        A_t(k,2*single_rob_len*nRobots+1+(k-1)*nRequired:2*single_rob_len*nRobots+k*nRequired) = c1*ones(1,nRequired);
%         
%        A_t(k,(k-1)*single_rob_len+1:k*single_rob_len) = c2*ones(1,single_rob_len);
%         
%        A_t(k,end) = -1;
%         
%     end
% 
%     b_t = zeros(nRobots,1);
%     
%     
%     A_t(nRobots+1,sol_len) = 1;
%     
%     b_t(nRobots+1) = 100;
    
    

end


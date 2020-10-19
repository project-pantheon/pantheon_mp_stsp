function T = calculateFinalTime(sol_edges,sel_vertex,Vr,c_req)

    %sel_vertex = sel_vertex(2:end); %take out depot

    t1 = 5*(sum(sol_edges)-2); %griglia sempre a 5x5
    
    %find stops
    [~,pos] = intersect(Vr,sel_vertex) ;
    
    t2 = sum(c_req(pos)); %peso sempre a 10
    
    
    T = t2 + t1;




end
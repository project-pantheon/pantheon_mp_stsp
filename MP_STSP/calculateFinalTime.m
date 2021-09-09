function T = calculateFinalTime(tour,sel_vertex,Vr,c_req)

    %sel_vertex = sel_vertex(2:end); %take out depot

    %t1 = 5*(sum(sol_edges)-2); %griglia sempre a 5x5
    t1 = 5*(size(tour,1)-3);
    
    %find stops
    [~,pos] = intersect(Vr,sel_vertex) ;
    
    t2 = sum(c_req(pos)); %peso sempre a 10
    
    
    T = t2 + t1;




end

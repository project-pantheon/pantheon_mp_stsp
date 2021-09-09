function tour = findTourReal(sel_edges,tot_required,sir)

%     sel_edges = round(edges(sol >= .99,:)); %necessary. Some variables are not 1

%    d_edges = findDoubledEdges(sel_edges);
    

%% TODO
%     tour = zeros(1,length(sel_edges));
%     
%     next_vertex = sel_edges(1,1);
% 
%     tour(1) = sel_edges(1,1); %always 1
%      
%     for i=2:length(sel_edges)+1 %include also the 1 at the end
%         
%        indx = find( sel_edges(:,1) == next_vertex);
%        tour(i) = sel_edges(indx,2);
%        next_vertex = sel_edges(indx,2);
%                        
%         
%     end

    % TOUR
    tour = zeros(1,length(sel_edges));
    mask = zeros(1,length(sel_edges));

    un = unique(sel_edges);
    
    for i = 1:length(un)
        cardinalita(i,:)= [un(i) length(find(sel_edges == un(i)))];
    end
    
    [tour, mask] = findTourRicReal(mask,sel_edges,cardinalita);
    
%     %first step
%     current = sel_edges(1,1);
%     tour(1) = current;
%     mask(1) = 1;
%     
%     i=2
%     while sum(mask)<=length(sel_edges)
%         new_elem = false
%         current = find(sel_edges(:,1)==sel_edges(current,2))
%         j=1;
%         while j<=size(current,1) || ~new_elem
%             if mask(current(j)) == 0 || ~new_elem
%                 tour(i) = sel_edges(current(j),1)
%                 current = find(sel_edges(:,1)==sel_edges(current,2))
%                 mask(i) = 1
%                 i=i+1
%                 new_elem=true
%             end
%         end
%     end

    tour = tour';
    
    tour(:,2) = zeros(length(tour),1);
    
    for i=sir
        count = 0;
        
        for j=1:length(tour)
           
        
            if (tour(j,1) == i) && (count == 0)
               
                tour(j,2) = 1;
                count = count + 1;
                
            end
        
        end
        
    end

    


end


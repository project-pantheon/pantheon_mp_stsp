function [tour,mask] = findTourRic(mask,sel_edges,cardinalita)
%FINDTOURRIC

%init
list_current = find(mask==0);
cnt=1;

origin_selected = false;

if isempty(mask)
    disp('A robot did not leave the depot!')
   return 
end

%origin selected using mask a reference
if sum(mask) == 0 %activated only one time
    origin=list_current(1);
else
    filter_edges=sel_edges.*mask';
    iter=1;
    while iter<=size(list_current,2) && ~origin_selected
        origin_in_prev_tour=find(filter_edges(:,1)==sel_edges(list_current(iter),1));
        if ~isempty(origin_in_prev_tour)
        	origin=list_current(iter);
            list_current= [list_current(iter) list_current(1:iter-1) list_current(iter+1:end)];
            origin_selected = true;
        end
        iter=iter+1;
    end
end

current=origin;

tour_closed = false;

% fill tour using mask
for i=1:size(list_current,2)
    if ~tour_closed
        tour(cnt) = sel_edges(current,1);
        cnt=cnt+1;
        mask(current)=1;
        current_candidates = find(sel_edges(:,1)==sel_edges(current,2))';
        
%         % Probably NOT required start____________________________________________________________________
%         %generate current_candidates_table to sort by cardinalita
%         current_candidates_table = []
%         for j=1:size(current_candidates,2)
%             current_candidates_table= [current_candidates_table; current_candidates(j), cardinalita(find(cardinalita(:,1)==sel_edges(current_candidates(j),2)),2)]
%         end
%         %sort current_candidates by cardinalita
%         current_candidates_table = sortrows(current_candidates_table,2)
%         %update current_candidates
%         current_candidates = current_candidates_table(:,1)'
%         % Probably NOT required end____________________________________________________________________
        
        
        
        candidate_selected=false;
        for j=current_candidates  % select candidate
            if mask(j)==0 && ~candidate_selected 
                current = j;
                candidate_selected=true;
            elseif j==origin %close tour
                tour(cnt) = sel_edges(origin,1);
                tour_closed = true;
            end
        end
    end
end

%merge tours
if sum(mask)<size(mask,2)
    [tour_,mask] = findTourRic(mask,sel_edges, cardinalita);
    
    if ~isempty(find(tour_==tour(1))) %tour_ contains the origin of tour
        tour = [tour_(1:find(tour_==tour(1))-1) tour tour_(find(tour_==tour(1))+1:end)];
    elseif ~isempty(find(tour==tour_(1))) %tour contains the origin of tour_
        tour = [tour(1:find(tour==tour_(1))-1) tour_ tour(find(tour==tour_(1))+1:end)];
    end
end
    

end


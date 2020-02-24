function [edges_sol] = find_edges(solution,edges)

edges_sol = edges(find(solution(1:length(solution))>0.1),:);

end


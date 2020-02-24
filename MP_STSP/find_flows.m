function [flows_sol] = find_flows(solution,edges)

flows_sol = [edges(find(solution > 0.1),:) solution(find(solution > 0.1))];



end


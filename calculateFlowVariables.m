function [delta_plus , delta_minus, Fa] = calculateFlowVariables(nStops,edges_list)

    %delta sets
    delta_plus = cell(1,nStops);
    delta_minus = cell(1,nStops);

    for ii = 1:nStops

        delta_plus{ii} = (edges_list(:,1) == ii);
        delta_minus{ii} = (edges_list(:,2) == ii);

    end

    % flow variables
    Fa = ones(length(edges_list),1)/((nStops)+1);


end


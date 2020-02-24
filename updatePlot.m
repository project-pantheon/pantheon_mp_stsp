function fig = updatePlot(robot,x_sol,edges,stops,color,Vri)

X = stops(:,1);
Y = stops(:,2);

segments = find(round(x_sol)); 

pt1 = zeros(3*length(segments),1);
pt2 = zeros(3*length(segments),1);

for ii = 1:length(segments)
    
    start = edges(segments(ii),1);
    stop = edges(segments(ii),2);
    
    pt1(3*ii-2:3*ii) = [X(start); X(stop); NaN];
    pt2(3*ii-2:3*ii) = [Y(start); Y(stop); NaN];  
end

figure(robot);

scatter(stops(:,1),stops(:,2),100,'o','filled','b')
scatter(stops(Vri,1),stops(Vri,2),100,'o','filled','r')
fig = plot(pt1,pt2,strcat('--',color),'LineWidth',2);

drawnow;
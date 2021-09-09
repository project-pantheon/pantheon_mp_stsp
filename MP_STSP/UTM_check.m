format long g

%Map
fnamemap = 'input/maps_config.json';
valmap = jsondecode(fileread(fnamemap));
map_index=0;

%Mission
fname = 'input/mission.json';
val = jsondecode(fileread(fname));

i=1;
while map_index==0 && i<=size(valmap,1)
    if strcmp(valmap(i).id,val.maps_config_id)
        map_index=i;
        id_mission_base=val.id_mission_base;
    end
    i=i+1;
end

%Trees
fname = "input/"+ valmap(map_index).mapping;
val_tab_trees = readtable(fname);
val_trees.id=val_tab_trees.Var1;
val_trees.name=[""];
i=1;
for i=1:size(val_tab_trees.Var2(:),1)
    val_trees.name(i)=val_tab_trees.Var2{i};
end

%Grid
fname = "input/"+ valmap(map_index).map_file;
%dlmwrite('filename.csv', [1.23456789], 'delimiter', ',', 'precision', 9)
val_grid = csvread(fname);
val_grid_utm=val_grid;
val_grid_utm_norm=val_grid_utm;


%UTM trasnformation
for i=1:size(val_grid_utm,1)
    [val_grid_utm(i,2),val_grid_utm(i,3),zone] = ll2utm(val_grid_utm(i,2),val_grid_utm(i,3));
    val_grid_utm_norm(i,2) = val_grid_utm(i,2)-val_grid_utm(1,2);
    val_grid_utm_norm(i,3) = val_grid_utm(i,3)-val_grid_utm(1,3);
end

scatter(val_grid_utm_norm(:,2),val_grid_utm_norm(:,3),10,'o','filled','b')

c = polyfit(val_grid_utm_norm(:,2),val_grid_utm_norm(:,3),1);

atan2(c(2),c1(1))
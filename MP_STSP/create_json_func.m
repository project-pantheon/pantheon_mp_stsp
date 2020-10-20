function [] = create_json_func(fname,nrobots,trees_rows,trees_cols,delta_rows,delta_cols,selected_id)
%CREATE_JSON requires all parameters of the experiment to create a json

fileID = fopen(fname, 'w+');

%build id struct
field = {   struct('trees_rows',num2str(trees_rows)), ...
            struct('trees_cols',num2str(trees_cols)), ...
            struct('delta_rows',delta_rows),  ...
            struct('delta_cols',delta_cols)};

%build id struct
trees = {};
for i_id=selected_id
    trees = {trees{:} struct('id', num2str(i_id)) };
end
full_trees = "0";

%struct build
struct_json.nrobots = num2str(nrobots);
struct_json.field = field;
struct_json.trees = trees;
struct_json.full_trees = full_trees;

%write json
text_json = jsonencode(struct_json);
fprintf(fileID,text_json);

fclose(fileID);

end


clear all
close all
clc

while 1
    
    % update data
    fname = 'exp_set/Set.json';
    set_struct = jsondecode(fileread(fname));
    fname = 'exp_set/Status.json';
    status_struct = jsondecode(fileread(fname));
    
    path_mission=sprintf("%s/M%d.json",set_struct.set{status_struct.set_status},status_struct.status);
    
    %change Mission file
    copyfile(path_mission,"Mission.json")
    
    %run algorithm
    MP_STSP_dev
    
    % update data
    fname = 'exp_set/Set.json';
    set_struct = jsondecode(fileread(fname));
    fname = 'exp_set/Status.json';
    status_struct = jsondecode(fileread(fname));
    
    path_data=sprintf("%s/M%d",set_struct.set{status_struct.set_status},status_struct.status);
    
    if solution_found
        save(path_data)
    end
    
    if status_struct.status+1 > status_struct.max
        status_struct.status = 1;
        if status_struct.set_status+1 > status_struct.set_max
            %reset set_status and break the loop
            status_struct.set_status = 1;
            fileID = fopen(fname, 'w+');
            text_status_json = jsonencode(status_struct);
            fprintf(fileID,text_status_json);
            fclose(fileID);
            break;
        else
            status_struct.set_status = status_struct.set_status+1;
        end
    else
        status_struct.status = status_struct.status +1;
    end
    
    
    fileID = fopen(fname, 'w+');
    text_status_json = jsonencode(status_struct);
    fprintf(fileID,text_status_json);
    fclose(fileID);
    
    
end

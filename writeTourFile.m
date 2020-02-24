function writeTourFile(tour)

%     sel_points = points(tour);
    fileID = fopen('/home/renzo/lab_ws/src/sherpa_ros/scripts/tour.csv','w');
    
    %fprintf(fileID,'id x y',tour');
    fprintf(fileID,'%d\n',tour');

    fclose(fileID);

end


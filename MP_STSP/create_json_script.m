clear all
close all

instances=10;

%nrobots = 3;
%trees_rows=8;
%trees_cols=8;
delta_rows=5;
delta_cols=5;

factors_str=["1_4","2_4","3_4","4_4"];
factors=    [ 1/4 , 2/4 , 3/4 , 4/4 ];

%factor_id = 1 ;
%1: 1/4
%2: 2/4
%3: 3/4
%4: 4/4


for nrobots=3:5
    
    for grid_size=4:2:14
        
        trees_rows=grid_size;
        trees_cols=grid_size;
        
        for factor_id=1:4
            
            files_directory= sprintf("exp_r%d_%dx%d_%s",nrobots,trees_rows+1,trees_cols+1,factors_str(factor_id));
            mkdir(files_directory)
            
            for i =1:instances
                
                fname = sprintf("%s/M%d.json",files_directory,i);
                
                %trees required
                n=round(trees_rows*trees_cols*factors(factor_id));
                
                % select randomly n unique ids
                rand_id = randi([1 trees_rows*trees_cols],1,n);
                rand_id = unique(rand_id);
                
                while (size(rand_id,2) < n)
                    rand_id = [rand_id randi([1 trees_rows*trees_cols],1,1)];
                    rand_id = unique(rand_id);
                end
                
                % selected_id
                selected_id = rand_id;
                
                create_json_func(fname,nrobots,trees_rows,trees_cols,delta_rows,delta_cols,selected_id)
            end
            
        end
        
    end
    
end
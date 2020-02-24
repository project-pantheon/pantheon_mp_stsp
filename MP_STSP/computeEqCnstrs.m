function [A, b] = computeEqCnstrs(idxs_complete, n_rows, n_cols)

n_turn_cnstrs =  (n_rows-2)*4+ (n_cols-2)*4 + (n_rows-2)*2*2  + (n_cols-2)*2*2 + 4;

n_idxs = length(idxs_complete);
A = sparse(zeros( n_turn_cnstrs, n_idxs*3 ));

debug_mat = zeros(n_rows, n_cols);
iter = 1;
                

for r = 1:n_rows
    r_inv = n_rows + 1 - r;
    for c = 1:n_cols
        
        if(r_inv==1)
            
            if(c==1)
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+1,n_cols*(r-1)+2]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+1, n_cols*(r-1) + 1 - n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 1;
            elseif(c==n_cols)
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+n_cols,n_cols*(r-1)+n_cols-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+n_cols, n_cols*(r-1)]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 1;
            else
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c+1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 2;
            end
            
        elseif(r_inv==n_rows)
            
            if(c==1)
                [id1, id1_rev] = findVector(idxs_complete, [1,2]);
                [id2, id2_rev] = findVector(idxs_complete, [1,n_cols + 1]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 1;
            elseif(c==n_cols)
                [id1, id1_rev] = findVector(idxs_complete, [n_cols,n_cols-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols,n_cols + n_cols]);
                %r, c, id1, id2_rev, id1_rev, id2
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 1;
            else
                [id1, id1_rev] = findVector(idxs_complete, [c,c - 1]);
                [id2, id2_rev] = findVector(idxs_complete, [c,c + n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [c,c + 1]);
                [id2, id2_rev] = findVector(idxs_complete, [c,c + n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 2;
            end
            
        else
            if(c==1)
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+1,n_cols*(r-1)+2]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+1,n_cols*(r-1)+1+n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+1,n_cols*(r-1)+2]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+1,n_cols*(r-1)+1-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 2;
            elseif(c==n_cols)
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r),n_cols*(r)-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r),n_cols*(r)+n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r),n_cols*(r)-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r),n_cols*(r)-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 2;
            else
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c+1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c+n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c+1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c+n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                
                [id1, id1_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-1]);
                [id2, id2_rev] = findVector(idxs_complete, [n_cols*(r-1)+c,n_cols*(r-1)+c-n_cols]);
                A(iter, id1) = 1; A(iter, id2_rev) = 1; A(iter, id1 + 2*n_idxs) = -1; iter = iter + 1;
                A(iter, id2) = 1; A(iter, id1_rev) = 1; A(iter, id2 + 2*n_idxs) = -1; iter = iter + 1;
                debug_mat(r_inv,c) = 4;
            end 
            
        end
    end
    
end

b = ones( iter-1, 1 );

end
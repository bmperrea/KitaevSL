function [var_con,var_dec,var_pos] = initseq_matlab_vec(mat)

    cols = size(mat,2);
    var_con = zeros(cols,1);
    var_dec = zeros(cols,1);
    var_pos = zeros(cols,1);
    
    for j = 1:cols
        [v_con,v_dec,v_pos] = initseq_matlab(mat(:,j));
        var_con(j) = v_con;
        var_dec(j) = v_dec;
        var_pos(j) = v_pos;
    end        

end
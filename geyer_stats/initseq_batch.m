function var_con = initseq_batch(x,nbatches)
    cols = size(x,2);
    rows = size(x,1);
    var_con = zeros(1,cols);
    batch_size = ceil(rows/nbatches);
    nbatches = floor(rows/batch_size);
    maxy = (nbatches*batch_size);
    pts = (1:maxy) + (rows-maxy); %prefer to sample the end of the chain
    for j=1:cols
        var_con(j) = initseq_vec(double(mean(reshape( x(pts,j), batch_size, nbatches),1)).');
    end
    var_con = var_con * batch_size;
end
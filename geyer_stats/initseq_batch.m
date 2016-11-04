function var_con = initseq_batch(x,nbatches)
    %Computes initseq more quickly by first reducing the data
    %The only output that is stored in this case in the best variance
    %estimator - the initial convex sequence estimator.
    %x is a matrix whose columns are Markov chains
    %nbatches is the number of batches. If the length of the chain is
    %smaller than nbatches the batch size is 1.
    cols = size(x,2);
    rows = size(x,1);
    var_con = zeros(1,cols);
    batch_size = ceil(rows/nbatches);
    nbatches = floor(rows/batch_size);
    maxy = (nbatches*batch_size);
    pts = (1:maxy) + (rows-maxy); %prefer to sample the end of the chain
    for j=1:cols
        var_con(j) = initseq_matlab(mean(reshape( x(pts,j), batch_size, nbatches),1).');
    end
    var_con = var_con * batch_size;
end
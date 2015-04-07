% computes the histogram for a laplacian of a grayscale image. 
function hist = get_texture_hist(bw_laplacian, bin_size)

    % init variables
    hist=[];
    min_t = 0;
    max_t = 255*8;
    
    nbins = ((max_t-min_t)/bin_size);
    bw_laplacian = bw_laplacian(:);
    
    % iteration over the bins
    for i=1:ceil(nbins)
        minbin=min_t+bin_size*(i-1);
        maxbin=min_t+bin_size*(i);
        
        % count the number of pizels that belong to current bin
        hist(i) = sum(bw_laplacian>=minbin & bw_laplacian<maxbin);
    end
end
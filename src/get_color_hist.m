function hist = get_color_hist(imgpix, bin_size)

    hist=[];
    r = imgpix(:,:,1);
    b = imgpix(:,:,2);
    g = imgpix(:,:,3);
    r = r(:);
    b = b(:);
    g = g(:);
    
    for i=1:ceil(255/bin_size)
        for j=1:ceil(255/bin_size)
            for k=1:ceil(255/bin_size)
                hist(i,j,k) = sum(r>=bin_size*(i-1) & r<bin_size*i... 
                    & b>=bin_size*(j-1) & b<bin_size*j... 
                    & g>=bin_size*(k-1) & g<bin_size*k);
            end
        end
    end
    
end
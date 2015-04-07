% converts a cell array of clusters to a grouping vector
function c_v= convert_cluster_cell_to_vec(clusts)
    c_v = zeros(1, 40);
    for clu=1:numel(clusts)
        c_v(clusts{clu}) = clu;
    end
end
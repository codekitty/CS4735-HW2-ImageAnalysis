% Performs a complete link or single link hierarchichal clustering down to
% n requested clusters. Returns C the final clustering and H the entire
% hierarchy.
function [C, H] = hier_clusters(A, n, type)

H = {};
D = 1-A; % convert from adjencency to distance
clusters = num2cell(1:size(A, 1)); % init clusters
while numel(clusters)>n
    nclusters = numel(clusters);
    I = zeros(nclusters, nclusters);
    
    % create cluster similarity matrix
    for i=1:nclusters
        for j=1:nclusters
            I(i,j) = compare_clusters(D, clusters{i}, clusters{j}, type);
        end
    end    

    % find the two closest clusters
    [I_sorted, iidx] = sort(I,2, 'ascend');
    [~, clostest_cluster] = min(sum(I_sorted(:,1:2), 2));
    if numel(clostest_cluster) > 1
        clostest_cluster = randsample(clostest_cluster,1);
    end
    
    % merge closest clusters
    merge_i = clostest_cluster;
    merge_j = iidx(merge_i, 2);
    clusters{merge_i} = union(clusters{merge_i},clusters{merge_j});
    clusters(merge_j) = []; % delete one cluster after merge
    H(end+1) = {clusters};
end

C = H{end};

end
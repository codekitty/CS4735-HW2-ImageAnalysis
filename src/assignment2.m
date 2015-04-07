% runs steps 1-4 of Assignment 2 in visual interfaces
function assignment2

close all;;
clc;

files = dir([curr_path filesep 'Images/*.ppm']);

imgObjs = cell(1,40);
bws = cell(1,40);
bws_laplacian = cell(1,40);
chists=[];
m_size = length(files); 
color_bin_size = 40;
texture_bin_size = 20;

% load and process all images by texture and color
for i=1:m_size
    imgObj = imread([curr_path filesep 'Images' filesep files(i).name]);
    imgObjs{i} = imgObj;
    
    % == ignore background and create color histograms
    r = imgObj(:,:,1);
    b = imgObj(:,:,2);
    g = imgObj(:,:,3);

    blackPixels = (double(r) < 25 & double(g) < 25 & double(b) < 25);
    [~,regions] = bwboundaries(blackPixels,'noholes');

    % sort by size
    table = tabulate(regions(:));   
    background_regions = table(table(:,2)>15,1);
    
    % remove continous regions of color
    rems = [];
    for regi=1:numel(background_regions)
        if any(~blackPixels(regions==background_regions(regi)))
            rems(end+1) = regi;
            break;
        end
    end
    background_regions(rems) = [];
    
    blackPixels = ismember(regions, background_regions);
    
    print_image_and_bin(imgObj, getInds(blackPixels), i, sprintf('%s%scolorbackgrounds25%s%g.png',curr_path, filesep, filesep, i));
    
    % set relevant pixels to nan
    r(blackPixels) = nan;
    b(blackPixels) = nan;
    g(blackPixels) = nan;  
    imgObj(:,:,1) = r;
    imgObj(:,:,2) = b;
    imgObj(:,:,3) = g;
    
    % get color histogram
    threedimhist = get_color_hist(imgObj, color_bin_size);
    chists(:,i) = threedimhist(:)./ sum(threedimhist(:));
    
    % == ignore background and create texture histograms
    bw = rgb2gray(imgObj);
    bws{i} = bw;

    % compute laplacian
    bw_laplacian = zeros(size(bw));
    for r=2:(size(bw(:,:), 1)-1)
        for c=2:(size(bw(:,:), 2)-1)
            bw_laplacian(r,c) =...
                double(bw(r,c))*9-sum(sum(double(bw(r-1:r+1, c-1:c+1))));   
        end
    end
    
    frac = 10;
    
    % remove background - if low texture and not in the center of image
    edgewidth  = floor(size(bw_laplacian,2)/frac);
    edgeheight = floor(size(bw_laplacian,1)/frac);

    % create dead zone mask
    mask = ones(1, size(bw_laplacian,2));
    mask(edgewidth:edgewidth*(frac-1)) = 0;

    % iterate over each image row
    for r=1:size(bw_laplacian,1)
        if (r>edgeheight && r<edgeheight*(frac-1))
            bw_laplacian(r, mask & abs(bw_laplacian(r, :)) < 3 ) = nan;
        else
            bw_laplacian(r, abs(bw_laplacian(r, :)) < 3 ) = nan;
        end
    end
    
	bws_laplacian{i} = bw_laplacian;
    
    print_image_and_bin(imgObj, getInds(isnan(bw_laplacian)), i, sprintf('%s%stexturebackgrounds25%s%g.png',curr_path, filesep, filesep, i));

    % compute texture histogram
    thists(:,i) = get_texture_hist(abs(bws_laplacian{i}),texture_bin_size);
    thists(:,i) = thists(:,i) ./ sum(sum(thists(:,i)));
end

% === Step 1: Gross color matching ==
C=compute_print_similarities(imgObjs, chists, 'color');

% === Step 2: Gross texture matching. ==
T=compute_print_similarities(imgObjs, thists, 'texture');


% === Step 3: Combine similarities, and cluster.
r=.6;
S = r*T + (1-r)*C;
print_similarities(imgObjs, S, sprintf('combination%gprctext',100*r));

% = Cluster to 7 clusters
complete_cluster = hier_clusters(S, 7, 'complete');
print_cluster_result(imgObjs,complete_cluster,'complete');


single_cluster   = hier_clusters(S, 7, 'single');
print_cluster_result(imgObjs,single_cluster,'single');

% === Part 4 benchmark against live human beings

[C1, T1, clusts1] = loadAndrey;
[C2, T2, clusts2] = loadRon;
[C3, T3, clusts3] = loadMichal;

print_img_similarity_mat(imgObjs, C, 'Friend1Color');
print_img_similarity_mat(imgObjs, T, 'Friend1Texture');
print_cluster_result(imgObjs,clusts,'Friend1');

print_img_similarity_idx(imgObjs, C2, 'Friend2Color');
print_img_similarity_idx(imgObjs, T2, 'Friend2Texture');
print_cluster_result(imgObjs,clusts2,'Friend2');

print_img_similarity_idx(imgObjs, C3, 'Friend3Color');
print_img_similarity_idx(imgObjs, T3, 'Friend3Texture');
print_cluster_result(imgObjs,clusts3,'Friend3');

% compare color and texture result to three friends
[~, idx] = sort(C, 2, 'descend');
CIDX = idx(:, [1:4 m_size-(2:-1:0)]);

ec_adj1 = compare_matches(CIDX, C1);
ec_adj2 = compare_matches(CIDX, C2);
ec_adj3 = compare_matches(CIDX, C3);

[~, idx] = sort(T, 2, 'descend');
TIDX = idx(:, [1:4 m_size-(2:-1:0)]);

et_adj1 = compare_matches(TIDX, T1);
et_adj2 = compare_matches(TIDX, T2);
et_adj3 = compare_matches(TIDX, T3);

fprintf('\nC vs. C friend1: %1.3f', ec_adj1);
fprintf('\nC vs. C friend2: %1.3f', ec_adj2);
fprintf('\nC vs. C friend3: %1.3f', ec_adj3);
fprintf('\nColor Similarity mean: %1.3f\n', mean([ec_adj1 ec_adj2 ec_adj3]));

fprintf('\nT vs. T friend1: %1.3f', et_adj1);
fprintf('\nT vs. T friend2: %1.3f', et_adj2);
fprintf('\nT vs. T friend3: %1.3f', et_adj3);
fprintf('\nTexture Similarity mean: %1.3f\n', mean([et_adj1 et_adj2 et_adj3]));

% compare cluster result to three friends 
c_v1 = convert_cluster_cell_to_vec(clusts1);
c_v2 = convert_cluster_cell_to_vec(clusts2);
c_v3 = convert_cluster_cell_to_vec(clusts3);

c_complete = convert_cluster_cell_to_vec(complete_cluster);
ec1 = clusteringError(c_v1, c_complete);
ec2 = clusteringError(c_v2, c_complete);
ec3 = clusteringError(c_v3, c_complete);

c_single = convert_cluster_cell_to_vec(single_cluster);
es1 = clusteringError(c_v1, c_single);
es2 = clusteringError(c_v2, c_single);
es3 = clusteringError(c_v3, c_single);

fprintf('\nR and Index - Complete link vs friend1: %1.3f', ec1);
fprintf('\nR and Index - Complete link vs friend2: %1.3f', ec2);
fprintf('\nR and Index - Complete link vs friend3: %1.3f', ec3);
fprintf('\nR and Index - Complete link mean: %1.3f\n', mean([ec1 ec2 ec3]));

fprintf('\nR and Index - Single link vs friend1: %1.3f', es1);
fprintf('\nR and Index - Single link vs friend2: %1.3f', es2);
fprintf('\nR and Index - Single link vs friend3: %1.3f', es3);
fprintf('\nR and Index - Single link mean: %1.3f\n', mean([es1 es2 es3]));
end

function print_cluster_result(imgObjs,clusters, folderprefix)

k=8; %number of pics to print in a row
for c=1:numel(clusters);
    cluster = clusters{c};
    
    print_all_imgs(imgObjs(cluster), cluster,...
      sprintf('%s-clusters/cluster%g.png',folderprefix,c),...
      k);
end

end

function print_img_row(imgObjs, inds,imgname, opt_h)

% optionally create a new figure
if ~exist('opt_h', 'var')
    m_size = numel(imgObjs);
    opt_h = figure('Position', [100, 100, 100*m_size, 80]);
end

% size of picture row
m_size = numel(imgObjs);

% print a row of images with a white label
for row=1:numel(imgObjs)
    subplot('position', [(row-1)/m_size, 0, 1/m_size, 1]);
    imshow(imgObjs{row});
    text(15,5,num2str(inds(row)), 'background', 'w');
end
fprintf('\n...Printing image to file %s', imgname);

% save to picture to file
screen2png(opt_h,imgname);

end

function print_all_imgs(imgObjs,inds,imgname, columns,opt_h)

rows = max(numel(imgObjs)/columns, 1);

% optionally create a new figure
if ~exist('opt_h', 'var')   
    opt_h = figure('Position', [100, 100, 100*columns, 80*rows]);
end

% print a row of images with a white label
for row=1:rows
    for column=1:columns
        if ((row-1)*columns+column > numel(imgObjs))
            break;
        end
        subplot('position', [(column-1)/columns, 1-((row)/rows), 1/columns, 1/rows]);
        imshow(imgObjs{(row-1)*columns+column});
        text(15,5,num2str(inds((row-1)*columns+column)), 'background', 'w');
    end
end
fprintf('\n...Printing image to file %s', imgname);

% save to picture to file
screen2png(opt_h,imgname);
end

function print_img_similarity_mat(imgObjs, S, foldername)
[~, idx] = sort(S, 2, 'descend');
SIDX = idx(:, [1:4 m_size-(2:-1:0)]);

print_img_similarity_idx(imgObjs, SIDX, foldername);
end

function print_img_similarity_idx(imgObjs, idx, foldername)

m_size = numel(imgObjs);
h = figure('Position', [100, 100, 700, 640]);
for row=1:8:m_size
    inds = idx(row:(row+8-1), :)';
    imgObj = imgObjs(inds(:));

    print_all_imgs(imgObj,...
                  inds(:),...
                  sprintf('%s%sImg%g-%g_similarity.png', foldername,filesep,row, row+8-1),...
                  7,...
                  h);
end

end

function print_image_and_bin(img,binimg,label, filename)
    h= figure(1);
    imshow(img);
    hold on;
    scatter(binimg(:,2), binimg(:,1), 40, '.b');
    hold off;
    text(15,5,num2str(label), 'background', 'w');

    
    % save to picture to file
    screen2png(h,filename);
end

function S=compute_print_similarities(imgObjs, hists, prefix)
S = pairwise_compare_hists(hists);

% print_similarities(imgObjs, S, prefix)
end

function print_similarities(imgObjs, S, prefix)
[S_sorted, idx(:,:)] = sort(S, 2, 'descend');

[~, l_sim] = min(sum(S_sorted(:,1:4),2));
[~, m_sim] = max(sum(S_sorted(:,1:4),2));

fprintf('\nMost similar group: %g', m_sim);

% most and least similar
imgname = sprintf('%s%smostleastsimilar/%smost.png',curr_path, filesep , prefix);
print_img_row(imgObjs(idx(m_sim,1:4)), idx(m_sim,1:4), imgname);
imgname = sprintf('%s%smostleastsimilar/%sleast.png',curr_path, filesep , prefix);
print_img_row(imgObjs(idx(l_sim,1:4)), idx(l_sim,1:4), imgname);

% print all color based similarities
output_folder = sprintf('%s%ssimilarities-%s',curr_path, filesep,prefix);
print_img_similarity_idx(imgObjs, idx(:, [1:4 40-(2:-1:0)]), output_folder);
end
% returns a two-column matrix with indices corresponding to the indices
% where the given matrix is positive
function inds=getInds(mat)
    [i,j,~] = find(mat);
    inds = [i,j];
end

function p = curr_path

    persistent f_p;

    if isempty(f_p)
        filename = mfilename('fullpath');
        [f_p, ~] = fileparts(filename);
    end
    p = f_p;
end
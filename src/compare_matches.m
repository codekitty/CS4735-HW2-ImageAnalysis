% compares two  matrices of that contain 3 similar and 3 unsimilar images
% per row. First column is ignored (self similarity) 
function accuracy = compare_matches(S, S1)

accuracy=0;
for i=1:size(S,1)

    % count a point if they share at least one similarity or dissimlarity
    if sum(ismember(S(i, 2:end), S1(i, 2:end)))>0
        accuracy = accuracy+1;
    end
end
accuracy = accuracy/size(S,1);
end
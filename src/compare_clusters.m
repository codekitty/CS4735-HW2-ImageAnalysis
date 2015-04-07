% compare clusters given adjecsency matrix
function s=compare_clusters(A, indsi, indsj, type)

if strcmp(type,'complete')
    s = max(max(A(indsi, indsj)));
else
    s = min(min(A(indsi, indsj)));
end

end
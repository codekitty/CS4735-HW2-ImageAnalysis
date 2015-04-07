function A=pairwise_compare_hists(hists)

m_size = size(hists,2);
A = zeros(m_size, m_size);  % distance matrix
for i=1:m_size
    for j=1:m_size
        A(i,j) = 1-l1_compare(hists(:,i), hists(:,j));
    end
end
end
% sum over the absolut difference at each vairable
function score = l1_compare(imghist1, imghist2)
    vec1 = imghist1(:);
    vec2 = imghist2(:);
    score = sum(abs(vec1./sum(vec1) - vec2./sum(vec2)))./2;
end
function success = cell2Histogram(cell)
    
% cell2Histogram with create a histogram of the data in a numerical 2D cell
% array made up of row vectors.

histData = [];

for i = 1:size(cell,1)
    for j = 1:size(cell,2)
        histData = [histData cell{i,j}];
    end
end

figure

histogram(histData);
title('Rays that hit Detector')
xlabel('Frequency')
ylabel('# of rays that hit')

success = 1;

end
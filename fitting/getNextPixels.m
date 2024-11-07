function nearestPixels = getNextPixels(processed, mask, fetchLen)
% getNextPixels - Finds the unprocessed pixels nearest to processed pixels
% 
%   nearestPixels = getNextPixels(processed, mask, fetchLen)
%
% Author:    Abbas El Moussawi
%
% DESCRIPTION:
% This function searches for the next target pixels to fit. It considers 
% two binary matrices represent-ing processed and masked pixels to identify 
% neighboring pixels within the specified region of inter-est. It begins by 
% isolating unprocessed pixels within the mask and those adjacent to processed
% ones using convolution with a Laplacian operator of a Gaussian function 
% to detect the edge. If adja-cent pixels are found, it calculates their nearest
% neighbors based on Euclidean distance. Otherwise, it locates the nearest
% unprocessed pixel within the mask. The output provides coordinates of the 
% nearest pixel, for extracting the fit parameters.
% the fit parameters. 
%
% INPUTS: 
%   processed      - A logic matrix indicating the processed pixels as true.
%   mask           - A logix matrix indicating pixels to ignore as false.
%   fetchLen       - The number of pixels to return.
%
% OUTPUT:
%   nearestPixels   - The nearest pixels found.
%
% SEE ALSO:
%  fitExcImageDatase, getNearestPixel, distEuc
% 
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen  


nearestPixels = [];
toFit = all(cat(3, ~processed,mask),3);

kernel = [0	0	1	0	0
0	1	1	1	0
1	1	0	1	1
0	1	1	1	0
0	0	1	0	0];
result = conv2(processed, kernel, 'same');
[y, x] = find(toFit == 1 & result > 0, fetchLen);

if ~isempty(x)
    nearestPixels = zeros(height(x),5);
    for i = 1:length(x)
        [nearest,d  ] = getNearestPixel(processed, y(i), x(i));
        nearestPixels(i,:)= [y(i) x(i) nearest d];
    end
elseif any(toFit(:))

    result = conv2(processed, kernel, 'same');
    [yy, xx] = find(toFit == 0 & result > 0 & result < 3, ceil(numel(toFit)/2));
    [tfy, tfx] = find(toFit == true);
    minD = Inf;
    minI = 1;
    for i = 1:height(yy)
        d = distEuc(tfy,tfx, yy(i), xx(i));
        md = min(d);
        if md < minD
            minI = i;
            minD = md;
        end
    end
    yy = yy(minI);
    xx = xx(minI);

    [nearest,d  ] = getNearestPixel(toFit, yy, xx);
    nearestPixels = [nearest yy xx d]; 
end
end

function [nearest, dist] = getNearestPixel(cond, refy,refx)


[y, x] = find(cond == 1);


dist = distEuc(refy,refx,y,x);

nearest = [y x dist];
nearest = sortrows(nearest,3);
nearest = nearest(~all(nearest == 0, 2), :);
dist = nearest(1, 3);
nearest = nearest(1,1:2);

end

function d = distEuc(y1,x1,y2,x2)
d = sqrt((x1 - x2).^2 + (y1 - y2).^2);
end

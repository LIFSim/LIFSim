function [imgs, n] = imagesLoader(folder,ext)
% imagesLoader - Retrieves images in a given folder sorted.
%
%   [imgs, n] = imagesLoader(folder,ext)
%
% Author: Abbas El Moussawi
%
% DESCRIPTION:
%   This loads images of a given extension with imread according to
%   ascending order of the name. E.g. 1.tif, 2.tif, ... 9.tif, 10.tif will
%   be loaded in some operating systems as 1.tif, 10.tif, 2.tif, ... This
%   function corrects this.
%
% INPUTS:
%   folder: The directory path of the folder where the images are stored.
%   ext:    The extension of the images, e.g. 'tif', 'tiff', etc..  
%
%
% OUTPUT:
%   imgs: And nx1 cell containing the images.
%   n: Number of images found and loaded
%   
%
% SEE ALSO:
%   loadTIFImages.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen
 
files= dir(fullfile(folder, ['*.' ext]));

n = length(files);
 
fileNumbers = zeros(n, 1);
for i = 1:n
    numStr = regexp(files(i).name, '\d+', 'match');
    fileNumbers(i) = str2double(numStr{1});
end

[~, sortedIndices] = sort(fileNumbers);
sorted = files(sortedIndices);
imgs = cell(n,1);


for i = 1:n
    imgs{i} = imread(fullfile(sorted(i).folder,sorted(i).name));
end

end


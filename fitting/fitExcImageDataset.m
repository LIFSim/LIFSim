function [fittedImage, fitResidue, tempFig] = fitExcImageDataset(...
    wnum, images3d, MM, linelist, startPos, startFit, mask,...
    options, plotImage, bounds, specParams, collParam)
% fitExcImageDataset - 
%
%   [fittedImage, fitResidue, tempFig] = fitExcImageDataset(...
%
% Author: Abbas El Moussawi and Torsten Endres
%
% DESCRIPTION:
%   Spatially resolved excitation spectra, named as excitation image dataset, 
% can be handled by this equation. The pixel positions are visited, and the 
% spectra are extracted from each corresponding position. To next positions to 
% handle are determined based on the previously fitted positions, here the function 
% getNextPixels is utilized. This function relies on an initial position to 
% start extract-ing the spectra from the image set and fitting. This initial 
% position should have been already fitted and the fit parameters must be provided. 
% Furthermore, a mask can be passed to this function to omit regions where no fit
% is required. This could be regions with low LIF intensity and regions out-side 
% the laser sheet. The number of simultaneous spectra to be fitted is determined by 
% the num-ber of available parallel workers. To use this function, configure the
% MATLAB parallel profile, as needed, by clicking the icon in the bottom left corner
% and then “Parallel Preferences”.
%
% INPUTS:
%   wnum: Wavenumber region in cm–1.  
%   images3d: The LIF image set stored in a 3D matrix with dimensions as follows:
%       -	x: Width of the image. 
%       -	y: Height of the image.
%       -	z: The image index. The length along this coordinate should be equal to 
%       wnum. For an im-age, idx, the image is then images3d(:,:,idx) that corresponds 
%       to the excitation wavenumber wnum(idx).
%   MM: Molar mass in g/mol.
%   linelist: The linelist obtained from function selectLines.
%   startPos: The (y, x) coordinates of the starting position in the image set. The 
%   spectrum at this position must be first fitted and the resulting parameters must be 
%   provided in startFit.
%   startFit: A matrix vector with seven places, which are the resulting fit parameters 
%   for startPos.
%   mask: A logical mask with height and width matching those of image3d. The positions 
%   set to true in this matrix would be processed.
%   options: The options objects needed by the optimization function lsqnonlin. Some 
%   presets can be obtained from the function getFitOptions.
% 
% OPTIONAL INPUTS:
%   plotImage: A logical variable to show the resulting temperature images as the fit
%   propagates through the positions (default = true). 
%   bounds: This defines the fit boundaries for the parameters. The boundaries here are
%   not in percentage. The upper and lower boundaries are then set around the result fit 
%   parameters of previously nearest fit. For the initial fit they are calculated from the 
%   startFit vector.
%   varargin: Other arguments to this function are not processed and would be redirected to fitExcitationSpec.
%  
% OUTPUT:
%   fittedImage: The resulting fit parameters are stored here. The height and width of 
%   this 3D ma-trix corresponds to image3d, and the z direction has the seven resulting 
%   fit parameters for each position.
%   fitResidue: A 2D matrix matching the height and width of image3d, storing the resulting 
%   residual of the fits.
%   timeElapsed: The duration of the fit in seconds.
%   tempFig: The figure object in case plotImage was set to true.   
%
% SEE ALSO:
%   fitExcitationSpec, parfeval, getFitOptions, getNextPixels, fitExcSpec-trumImage.mlx, sensitivityAnalysis.mlx
%
% COPYRIGHT 2024:
%   EMPI-RF - University of Duisburg-Essen    

arguments
    wnum (1,:) double
    images3d (:,:,:) double
    MM (1,1) double
    linelist (:,:) cell
    startPos (1,2) double
    startFit (1,:) double
    mask (:,:) logical
    options (1,1)
    plotImage (1,1) logical = true
    bounds (1,:) double = [1000   1   1   0.15   1   1   0.07]
    specParams.resFactor (1,1) double = 3
    specParams.Z (:,:) double = [startFit(1) 1;startFit(1)+1 1]
    specParams.limit (1,1) double  = 12

    collParam.colls = {}
    collParam.gas = {}
    collParam.P (1,1) double = 1
end
  
starty = startPos(1);
startx = startPos(2);

[h, w, ~]= size(images3d);
options.Display="none";
fitted = false(h, w);
fitted(starty,startx) = true;
fitResidue = NaN(h,w);
fittedImage = Inf(h,w,length(startFit));


measSpec = images3d(starty,startx,:);
measSpec = measSpec(:);
measSpec = measSpec./max(measSpec);

lb = startFit-bounds;
ub = startFit+bounds;

p =  namedargs2cell(specParams);


c = namedargs2cell(collParam);

[startFit, res] = fitExcitationSpec(wnum, measSpec, linelist, MM, ...
    startFit,options,p{:}, c{:}, lb=lb, ub= ub);



fittedImage(starty, startx, : ) = startFit;
fitResidue(starty, startx, : ) = res;

pixelsToFit = length(mask(mask==true));


p = gcp();


simultaneous = p.NumWorkers;
tempFig = {};
if plotImage
    tempFig = figure;
    tempFig.WindowState = 'maximized';

end

for px = 1:pixelsToFit
    fitPool = getNextPixels(fitted, mask, simultaneous);
    if isempty(fitPool)
        break;
    end
    lenPool = size(fitPool,1);
    fEval(1:lenPool) = parallel.FevalFuture;
    for idx = 1:lenPool

        y = fitPool(idx, 1);
        x = fitPool(idx, 2);

        yf = fitPool(idx, 3);
        xf = fitPool(idx, 4);

        measSpec = images3d(y,x,:);
        measSpec = measSpec(:);
        measSpec = measSpec./max(measSpec);

        nearestFit = fittedImage(yf, xf, : );
        nearestFit = nearestFit(:);

        lb = nearestFit-bounds;
        ub = nearestFit+bounds;

        pxfitHandle = @(~) fitExcitationSpec(...
            wnum, measSpec, linelist, MM, nearestFit(:), options, ...
            resFactor=specParams.resFactor, Z=specParams.Z, limit=specParams.limit,...
            colls=collParam.colls, gas=collParam.gas, P=collParam.P,  lb=lb, ub= ub);

        fEval(idx) = parfeval(p, pxfitHandle, 2);


    end

    for i = 1:lenPool
        [completedIdx,value, res] = fetchNext(fEval);
        y = fitPool(completedIdx, 1);
        x = fitPool(completedIdx, 2);
        fittedImage(y, x, : ) = value;
        fitResidue(y,x) = res;
        fitted(y,x) = true;
    end
    if plotImage
        figure(tempFig);
        imagesc(fittedImage(:,:,1));
        axis image;

        title('Temperature image')
        colorbarJet(tempFig, label='Temperature / K');

        drawnow;

    end
end


if plotImage

    figure(tempFig);
    t = fittedImage(:,:,1);
    axis image;

    t(t==0)=nan;
    imagesc(t);
    axis image;
    title('Temperature image')
    colorbarJet(tempFig, label='Temperature / K');

    drawnow;
end

end

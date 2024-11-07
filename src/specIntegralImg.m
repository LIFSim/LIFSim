function img = specIntegralImg(fitImage, images3d, linesCell,...
                                 lineName, boundaries, scanStep,MM, varargin)
% specIntegralImg(fitImage, images3d, lines_use)
% Calculates pixelwise the spectrum integral and returns a 
% a corresponding image.
% 
%  Author:     Abbas El Moussawi | 29.01.2020
% 
%  Data for 'specIntegralImg':
%       fitImage        : The fit values from lifmat (known at T_image).
%       images3d        : The lifmat input binned images in one 3d array.
%       lines_use       : Lines table used in lifmat.
%       lineName        : Which peak to simulate, e.g. input str 'Q(42)'
%       boundaries      : How far in wavenumbers to simulate around the peak 
%                         e.g. Q(42)-1 till Q(42)+4. E.g. input [-1 4]
%       scanStep        : Steps -1:scanStep:4. E.g. input 0.025
%  Output:
%       img             : img = pixelwise integeral of spectrum

    for y=1:(size(images3d,1))
        for x=1:(size(images3d,2))
            scalingimg(y,x)= fitImage(y,x,4).*max(images3d(y,x,:));
        end
    end

    
    line = getLineInfo(linesCell, lineName);

    wnum=((line{4}+boundaries(1)):scanStep:(line{4}+boundaries(2)))';

    img = zeros(size(images3d(:,:,1)));
    lines = cell2mat(linesCell(:,2:10));

    for y=1:(size(images3d,1))
        for x=1:(size(images3d,2))
            fitVals = squeeze(fitImage(y,x,:));
            spec = excitationSpec(...
                wnum,...
                lines,...
                MM, ...
                fitVals(1),...%T
                fitVals(2),... %dnuGL
                fitVals(3),... %dnuLL
                'dnuL',fitVals(4),... %dnuL
                'scale',fitVals(5),... %scale
                'dnuSh',fitVals(6),... %shift
                'offset',fitVals(7),... %offset
                varargin{:}...
                );
%             spec = scalingimg(y,x).*lifmat_spectrum(...
%                 wnum,...
%                 linesCell,...
%                 fitImage(y,x,1),...
%                 fitImage(y,x,2),...
%                 fitImage(y,x,3),...
%                 fitImage(y,x,4),...
%                 fitImage(y,x,5),...
%                 fitImage(y,x,6)...
%                 );
            img(y,x) = sum(spec(:).*scanStep); 
        end
    end

end
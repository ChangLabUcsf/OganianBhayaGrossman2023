function [cm2d] = ecog_2dcm(nColors, edgeC)
% create a 2 dimensional colormap with four different edge colors and white
% in center.
if nargin<1 
    nColors = 100; % colormap resolution
end

if nargin < 2 % edge colors. middle is always white
    edgeC = [0 0 1; 1 0 0; 200/256 0 200/256; 0 200/256 0];
end
% single colors
cmaps = cell(4,1);
for i = 1:4
    cmaps{i} = zeros(nColors,3);
    for cc = 1:3
        cmaps{i}(:,cc) = linspace(edgeC(i,cc),1, nColors)';
    end
end
% combine to colormap
fullmap = cell(2,1);
for i = 1:2
    fullmap{i} = [cmaps{(i-1)*2+1}; cmaps{i*2}(end:-1:1,:)];
end

% combine into 2d map
cm2d = zeros(nColors*2,nColors*2,3);
for i = 1:nColors*2
    for j = 1:nColors*2
        cm2d(i,j,:) =  fullmap{1}(j,:)*0.5 +  fullmap{2}(i,:)*0.5;
    end
end

if nargin<3
    figure;
    imagesc('XData', -nColors:nColors, 'YData',  -nColors:nColors, 'Cdata', cm2d); 
    axis tight;
end

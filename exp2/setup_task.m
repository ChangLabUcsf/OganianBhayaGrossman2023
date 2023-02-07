clc; clear all;

stimFolder = '../../data/exp2/stim_info';
outpath = '../data/exp2';
addpath('util');


%% plot defitions
nColors = 100;
cols = 'rbgmc';
rbcm = [[repmat(linspace(0,1, nColors)',1,2),ones(nColors, 1)] ; [ones(nColors, 1), repmat(linspace(1,0, nColors)',1,2)]];
rcm = [ones(nColors, 1), repmat(linspace(1,0, nColors)',1,2)];
cm2d = ecog_2dcm;
%% pink color map
pinkcurmap = pink(nColors); pinkcurmap = pinkcurmap(end:-1:1,:);
%% new colormap
rbcm = [[repmat(linspace(0,1, nColors)',1,2),ones(nColors, 1)] ; [ones(nColors, 1), repmat(linspace(1,0, nColors)',1,2)]];
newc1 = [1 1 1];%[202 60 172]/256;
newc1 = [120 174 220]/256;
% newc2 = [245,245,220]/256;
newc2 = [1 1 1];%[0 0 0];
nColors=100;
newmap = [linspace(newc2(1), newc1(1), nColors)',linspace(newc2(2), newc1(2), nColors)',linspace(newc2(3), newc1(3), nColors)'];
newmap(1,:) = 1 -newc2;%[1 1 1];
colormap(newmap)
axis off
%% vowel space colors
vowSpCol = [0.5 0.5 0.5;0 0 0];

%% dimex contour
dimCont = load('dimexContour.mat');
figure;
dimCont.M = contour(dimCont.C{1}, dimCont.C{2},dimCont.N', [50 50]);
%% vowels 
vowNames = {'a', 'e', 'i', 'o', 'u'};

%% participant list
SID = {'EC214'; 'EC219'; 'EC221'; 'EC222'; 'EC225'; 'EC235'; 'EC242'; 'EC243'};

allHemi = {'lh', 'lh', 'lh', 'lh', 'lh' ,'lh', 'rh', 'lh'};

taskname = 'task_bilingVowelSpace';


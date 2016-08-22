clear all; close all; clc;

load('fhr.mat');
fhr = data;
M = 24; % last 24 minutes
T = 2; % time series length (in minutes)
fs = 4; % 4Hz
numMinsDelete = 5; % delete last 5 minutes
overlap = 1; % 1 minute of overlap, 1 minute of new data

featureList{1, 1} = 'MedianValue';
featureList{2, 1} = 'MedianDeviation';
featureList{3, 1} = 'stv';
featureList{3, 2} = fs;
featureList{4, 1} = 'sti';
featureList{4, 2} = fs;
featureList{5, 1} = 'msti';
featureList{5, 2} = fs;
featureList{6, 1} = 'poincare';
featureList{6, 2} = [1];
featureList{7, 1} = 'higuchi';
featureList{8, 1} = 'EstPsd';
f = [.06 .3 1];
featureList{8, 2} = [.06 .3 1 fs];
featureList{9, 1} = 'ltv';
featureList{10, 1} = 'delta';
featureList{10, 2} = [fs];
featureList{11, 1} = 'SampEn';
featureList{11, 2} = [2];
featureList{12, 1} = 'FuzzyEn';
featureList{12, 2} = [2, 2];
featureList{13, 1} = 'additionalFeatures';
featureList{13, 2} = fs;

computedFeatures = getFHRfeatures(fhr, featureList, M, T, overlap, fs, numMinsDelete);

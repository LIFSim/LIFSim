rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));
cd(rootDir);

clearvars -except rootDir
format LONG G
format COMPACT
close all
clc

availableSpecies = getAvailableSpecies;

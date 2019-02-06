%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released under: MIT License
% 2019, Sebastiano Vascon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('dynamics'));

fprintf('Compiling Infection-Immunization dynamics...\n');
mex dynamics/inImDynC.cpp -outdir dynamics
fprintf('done!\n');
    

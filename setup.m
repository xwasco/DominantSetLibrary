addpath(genpath('dynamics'));

if exist('inImDynC','file')~=3
    fprintf('Compiling Infection-Immunization dynamics...');
    mex dynamics/inImDynC.cpp -outdir dynamics
    fprintf('done!\n');
end
    
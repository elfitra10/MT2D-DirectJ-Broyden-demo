function run_full_inversion()
%% run_full_inversion.m
% Run the English clean release version with model, RMS, and data-fit figures.

clc; close all;

repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root)
    repo_root = pwd;
end
cd(repo_root);

addpath(genpath(fullfile(repo_root, 'src')));

main_script = fullfile(repo_root, 'src', 'Inversi_TETM_DirectJ_Broyden.m');
if exist(main_script, 'file') ~= 2
    error('Main English clean-plots inversion script not found: %s', main_script);
end

fprintf('\nRunning English clean MT2D DirectJ-Broyden demo with figures.\n');
fprintf('Repository root: %s\n', repo_root);
fprintf('Main script    : %s\n\n', main_script);

run(main_script);
end

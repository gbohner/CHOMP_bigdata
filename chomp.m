function [opt, ROI_mask, ROIs] = chomp( opt )
%CHOMP This function automatically extracts regions of interest from a 
% two-photon microscopy recording of neuronal activity visualized by calcium-reporters.
% 
%   opt is a class of chomp-options
%
%   Written by Gergo Bohner <gbohner@gatsby.ucl.ac.uk> - 2015/11/12
%   Last update - 2016/03/15

gtic = tic;

%Find directory we're running from and add all subfolders to matlab path
dir_orig = pwd;
opt.code_path = [fileparts(mfilename('fullpath')) filesep];
addpath(genpath(opt.code_path));
cd(opt.code_path);

%Set current timestamp if not provided
if isempty(opt.timestamp) %else use the specified timestamped file for re-analysis
  opt.timestamp = datestr(now, 30); 
end 

%Make sure to create folders required;
[s,mess,messid] = mkdir([opt.root_folder opt.input_folder]);
[s,mess,messid] = mkdir([opt.root_folder opt.output_folder]);
[s,mess,messid] = mkdir([opt.root_folder opt.precomputed_folder]);
[s,mess,messid] = mkdir([opt.root_folder opt.results_folder]);


%Stabilize the image stack
opt = StabilizeFrames(opt);
toc(gtic);

%Preprocess the data, store it in the input folder, time-stamped
opt = extractData(opt);
toc(gtic);

%Learn the cell model from the data, store in the output folder
Model_learn(opt);
toc(gtic);

%Collect results
[ROI_mask, ROIs] = getROIs(opt);
toc(gtic);

%Clean up if required (delete unnecessary inbetween files etc)
if opt.cleanup
end

cd(dir_orig);

end


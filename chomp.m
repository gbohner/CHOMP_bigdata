function [opt, ROI_mask, ROIs] = chomp( opt )
%CHOMP This function automatically extracts regions of interest from a 
% two-photon microscopy recording of neuronal activity visualized by calcium-reporters.
% 
%   opt is a class of chomp-options
%
%   Written by Gergo Bohner <gbohner@gatsby.ucl.ac.uk> - 2015/11/12

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


%Stabilize the data if asked for:
if opt.stabilize
  opt = StabilizeFrames(opt);
end

%Preprocess the data, store it in tmp folder, time-stamped. Make sure to
%return updated opt for futher processing
opt = extractData(opt);

%Learn the cell model from the data
Model_learn(opt);

%Collect results
[ROI_mask, ROIs] = getROIs(opt);

%Clean up if required
if opt.cleanup
end

cd(dir_orig);

end

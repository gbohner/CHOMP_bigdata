function [opt, ROI_mask, ROIs] = chomp( opt )
%CHOMP This function automatically extracts regions of interest from a 
% two-photon microscopy recording of neuronal activity visualized by calcium-reporters.
% 
%   opt is a struct containing field-value pairs of options, 
%   see the opt_default file for possible fields and explanations
%
%   Written by Gergo Bohner <gbohner@gatsby.ucl.ac.uk> - 2015/11/12

%Find directory we're running from and add all subfolders to matlab path
dir_orig = pwd;
opt.code_path = [fileparts(mfilename('fullpath')) filesep];
addpath(genpath(opt.code_path));
cd(opt.code_path);

%Get current timestamp, we'll refer to all intermediate files via this
opt_def = opt_default();
opt_def.timestamp = datestr(now, 30);

%Generate opt structure via loading in defaults, and overwriting with the ones coming from input
opt = struct_merge(opt_def, opt);

%Make sure to create folders required;
[s,mess,messid] = mkdir(opt.input_folder);
[s,mess,messid] = mkdir(opt.output_folder);
[s,mess,messid] = mkdir(opt.precomputed_folder);


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


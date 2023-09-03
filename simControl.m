function simControl = simControl()

global example;

example = 2;
% set to 1 or 2 to run the supplied examples. Used in setTarget.m and simControl.m
% Set to 0 when not using example files (user must modify simControl.m and setTarget.m)


%% user-supplied file paths
spicePath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe'; % This is the path to your LT Spice installation
filePath = 'C:\Users\radam\Documents\LTspiceXVII\'; %This is the path to the working LTSPICE folder (schems, netlists, simulation output files)
% ******************************************************************* %

%% Examples 1 or 2 (see readme) **************************************

fprintf('example = %d\n',example);
% ******************************************************************* %
if example==1 % 3rd-order 1 op-amp lowpass
    fileName = 'example1'; % name of the LTSpice schematic you want to optimize (without the .asc)
end
if example==2 % weird elliptic filter, 1 op-amp with an inductor
    fileName = 'example2'; % name of the LTSpice schematic you want to optimize (without the .asc)
end

% **** if example set to 0, user must fill in their own filename below ****
%   fileName = 'yourLTspiceFilename'; % note, leave off the .asc extension


%% Cell arrays that are filled out by the user;

%% simControlOPtInstNames - 
%      cell array containing the instance names of the components the
%      optimizer will control. Enclose each instance name in single quotes


%% simControlMinVals 
        % The min value for the above instances that the the optimizer is
        % allowed to use. Note, use standard numerical entry, avoid using
        % 'k' or 'p' or 'n', etc.

%% simControlMaxVals - 
        % Cell array containing the min value for the above instances that the the optimizer is allowed to use

%% simControlInstTol
        % Component tolerances, applied post-optimization (during the schematic generation process)
        % Values enetered in 'E' format as below (enclose each in single quotes)
        % Applies only to instances that are being optimized
 
        % 20%    E6
        % 10%    E12
        % 5%     E24
        % 2%     E48
        % 1%     E96
        % 0.5%   E192
        % detailed explanation; <https://en.wikipedia.org/wiki/E_series_of_preferred_numbers>



%% example circuits, see README. Use these as a reference for your own circuits

if example==1
    simControlOPtInstNames = {'R2'    'R4'   'C1'   'C2'   'C3'};
    simControlMinVals =      {100     100    1e-12  1e-12   1e-12};
    simControlMaxVals =      {1e5     2.5e3  1e-6   1e-6   1e-6};
    simControlInstTol =      {'E96'   'E96'  'E24'  'E24'  'E24'};
    LTSPice_output_node = 'V(vout)';
end

if example==2
    simControlOPtInstNames = {'R3'   'R2'    'R4'   'C1'   'C2'   'C3'  'L1'  'C4'};
    simControlMinVals =       {1000   100     1000   1e-12  1e-12  30e-12 1e-5  1e-12}; % note, use numeric values only (no k, u, etc)
    simControlMaxVals =       {10e3    1e5     2.5e3  1e-6   1e-6   1e-6  1e-3  1e-5}; % note, use numeric values only (no k, u, etc)
    simControlInstTol =       {'E96'    'E96'   'E96'  'E24'  'E24'  'E24'  'E12'  'E24'};
    LTSPice_output_node = 'V(vout)';
end





%% set the match mode. 
% 1 = amplitude only
% 2 = phase only
% 3 = amplitude and phase both

matchMode = 1;

% set return list
simControl = {fileName spicePath filePath simControlOPtInstNames simControlMinVals simControlMaxVals simControlInstTol LTSPice_output_node matchMode};

end
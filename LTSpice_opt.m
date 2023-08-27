clear all;
close all;
fclose('all');


global passCell; % used to pass vars to the opt function without needing lots of individual globals


example = 1; 
% set to 1 or 2 to run the supplied examples. Used in setTarget.m and simControl.m
% passed to simControl. 
% Set to 0 when not using example files (user must modify simControl.m and setTarget.m)


% sources we are borrowing from;
% https://medium.com/@amattmiller/running-ltspice-from-matlab-630d551032cc.
% https://www.mathworks.com/matlabcentral/fileexchange/48840-round-to-electronic-component-values



simControl = simControl(example); % read simulation control input file, filled out by user

fileName = simControl{1};
spicePath = simControl{2};
filePath = simControl{3};
simControlOPtInstName = simControl{4};
simControlMinVals =  simControl{5};
simControlMaxVals = simControl{6};
simControlInstTol = simControl{7};
simControlOPtInstNames = simControl{4};
LTSpice_simTime = simControl{8}; % how long to wait for LTSPice to finish sim
LTSpice_output_node = simControl{9};
matchMode = simControl{10}; % 1 = ampl only, 2 = phase on;y, 3=both

% the current Matlab files are in dropbox/matlab_cornell/matlab_ltspice

% derived file paths and run scripts
netlist_fname = sprintf('%s%s.net',filePath,fileName); % netlist filename
RunLTstring = sprintf('start "LTspice" "%s" -b "%s%s.net"',spicePath, filePath, fileName);
LTSpice_outputfile = sprintf('%s%s.raw', filePath,fileName);


passCell = cell(18,1); % used to pass lots of variables at once into functions
passCell{1} = spicePath;
passCell{2} = filePath;
passCell{4} = fileName;

passCell{10} = netlist_fname;
passCell{12} = RunLTstring;
passCell{13} = LTSpice_outputfile;
passCell{16} = LTSpice_simTime;
passCell{17} = LTSpice_output_node;
passCell{18} = matchMode;

% send command to write netlist
string = sprintf('"%s" -netlist "%s%s.asc"',spicePath, filePath, fileName);
fprintf('Issuing command to write LTspice netlist\n%s\n',string);

[status,cmdout] = system(string);

pause(LTSpice_simTime);



% run initial simulation once to get frequencies, so we can computed the
% target response on the same freqyency grid as the simulation
fprintf('Issuing command to run initial LTspice simulation\n%s\n',RunLTstring);

[status,cmdout] = system(RunLTstring);
pause(LTSpice_simTime);
if(status) 
    fprintf('ERROR, LTSpice sim failed to run. Check setup or increase LTSpice_simTime\n ');
    return;
end

try
    result = LTspice2Matlab(LTSpice_outputfile); % parse output file
catch
    fprintf('ltspice2matlab error try again in 1/2 sec!\n');
    pause(0.5);
    result = LTspice2Matlab(LTSpice_outputfile); % parse output file
end

freqx = result.freq_vect;
passCell{15} = freqx;

% plot initial response
found_node=0;
for i = 1:result.num_variables
    if strcmp(result.variable_name_list{i},LTSpice_output_node)
        found_node=1;
        fresp = abs(result.variable_mat(i,:));
        if(matchMode==2 || matchMode==3)
           phase = unwrap(angle(result.variable_mat(i,:)));
        end
        if matchMode==1 % ampl only match
            figure;
            semilogx(freqx,20*log10(fresp));
            title('initial LTspice ampl response (dB) before optimization');
        end
        if matchMode==2 % phase only match
            figure;
            semilogx(freqx,phase);
            title('initial LTspice phase response (radians) before optimization');
        end
        if matchMode==3 % both phase and ampl match
            figure;
            subplot(1,2,1),semilogx(freqx,20*log10(fresp)); 
            title('initial LTspice freq response (dB) before optimization');
            subplot(1,2,2),semilogx(freqx,phase); 
            title('initial LTspice phase response (radians) before optimization');
        end
    end
end
if(found_node==0)
    fprintf('ERROR, output node %s not found in netlist\n',LTSpice_output_node);
    return;
end

% set the target response, using the frequencies that match the LTSpice sim (freqx)

temp_return = setTarget(passCell,freqx,matchMode,example);
target = temp_return{1};
errWeights = temp_return{2};

passCell{5} = target;
passCell{3} = errWeights;

% read in initial netlist

fid = fopen(netlist_fname);
netlist = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
netlist = netlist{1};
netlist = regexp(netlist,' ','split'); % netlist(n) contains the entire line; and each line contains a cell array, with every entry in a different cell
numlines_netlist = size(netlist,1);

passCell{9} = netlist;
passCell{11} = numlines_netlist;

% find how many compnents are being Optd, and make an index that points
% to the line number of those components
numOptd = size(simControlOPtInstName,2); % number of instances that are being optimized
OptLine = zeros(numOptd,1); % an array that points to the netlist lines where the inst name contains 'x'
UB = zeros(numOptd,1); % upper bound for optimizer
LB = zeros(numOptd,1); % lower bound for optimizer

kkk=1;
for k = 1:numlines_netlist % go through lines
    N=size(netlist{k},2);
    thisLine = netlist{k};
    % for each netlist line, check all the instance names that are being
    % optimized, and if there is a match, save the line # in the netlist
    % file so the optimizer knows what to vary
    for kk = 1:numOptd
        if contains(thisLine{1},char(simControlOPtInstName(kk)))
            OptLine(kkk)=k;
            UB(kkk) = cell2mat(simControlMaxVals(kk)); % upper bound to pass to optimizer
            LB(kkk) = cell2mat(simControlMinVals(kk)); % upper bound to pass to optimizer
            kkk = kkk + 1;
         
        end
    end
end

numMatchingInstFound=kkk-1;

if(numOptd ~= numMatchingInstFound)
    fprintf('ERROR;\n');
    fprintf('number of instances to be optimized in control file = %d\n',numOptd);
    fprintf('number of matching instances found in netlist = %d\n',numMatchingInstFound);
    fprintf('check Instance name spelling in control file\n');
    return;
end

      
        

passCell{6} = numOptd;
passCell{7} = OptLine;



% create and initialize the optParams array which is passed to the optimizer
%optParams = zeros(numOptd,1); % this holds the current optimizer params, init to schematic values with 'xx' in inst name
nomParams = zeros(numOptd,1); % this holds the nominal values init to schematic values with 'xx' in inst name, doesnt change

for k=1:numOptd
    thisLine = netlist{OptLine(k)}; % only lines that will be optimized here
    newStr = char(thisLine{4}); % in case its just a number
    % replace any 'micro' symbols from LTSpice with 'u'
    for kk=1:length(newStr)
        temp = uint32(newStr(kk));
        if(temp==181) % micro symbol, looks like 'u' but is not!!
            temp=117; % replace with an actual 'u'
            newStr(kk) = temp;
        end
        
    end


    test = str2double(newStr);  %if this is nan then we must substitute for the k or u
    if isnan(test) % must contain k or K or M or pf or nf or uf or mf
       
        % newStr = char(thisLine{4});
        newStr = strrep(newStr,'M','e6');
        newStr = strrep(newStr,'G','e9');
        newStr = lower(newStr);
        newStr = strrep(newStr,'k','e3');
        newStr = strrep(newStr,'pf','e-12');
        newStr = strrep(newStr,'ph','e-12');
        newStr = strrep(newStr,'p','e-12');
        newStr = strrep(newStr,'nf','e-9');
        newStr = strrep(newStr,'nh','e-9');
        newStr = strrep(newStr,'n','e-9');
        newStr = strrep(newStr,'uf','e-6');
        newStr = strrep(newStr,'uh','e-6');
        newStr = strrep(newStr,'u','e-6');
        newStr = strrep(newStr,'mf','e-3');
        newStr = strrep(newStr,'mh','e-3');
        newStr = strrep(newStr,'m','e-3');


    end


    nomParams(k) = str2double(newStr); % assumes for R's and C's that the value is the 4th entry (inst, node, node, value)

end

passCell{8} = nomParams;

fprintf('\n****************\n**************\nEntering Optimization Loop, please be patient ...\n************\n***********\n');
UB = UB./nomParams; % translate upper bounds into relative upper bounds
LB = LB./nomParams; % translate lower bounds into relative lower bounds
optParams = ones(numOptd,1); % the params used by the optimizer are multiplicative factors applied to the original components (in nomParams)

X = lsqnonlin(@optLTspice,optParams, LB, UB); %************************************ OPT OPT OPT *******************

passCell{14} = X;

fprintf('\n*************\n************\nDONE! Generating outputs ...\n***********\n*********\n');

for k=1:numOptd
    fprintf('%s %2.12e\n',char(netlist{OptLine(k)}(1)),X(k)*nomParams(k));
   
end


% re-run sim with current netlist

system(RunLTstring);
pause(LTSpice_simTime); 

result = LTspice2Matlab(LTSpice_outputfile); % parse output file
freqx = result.freq_vect;
pause(0.1);

% get opt response, no quantization
for i = 1:result.num_variables
    if strcmp(result.variable_name_list{i},LTSpice_output_node)
        fresp_opt = abs(result.variable_mat(i,:));
        phase_opt = unwrap(angle(result.variable_mat(i,:)));
    end
end


update_schematic(passCell,simControl); % write a new _opt schematic. Quantization is applied on the write-out
pause(0.1);
fprintf('\n\nNew schematic with optimum component values generated\nFilename = %s%s_opt.asc\n\n',filePath, fileName);
% re-run sim on the "_opt" schematic to check the quantization

opt_schem_fname = sprintf('%s%s_opt.asc',filePath,fileName); % opt netlist filename
RunLTstring_opt = sprintf('start "LTspice" "%s" -b "%s%s_opt.net"',spicePath, filePath, fileName);
LTSpice_outputfile_opt = sprintf('%s%s_opt.raw', filePath,fileName);


% send command to write netlist
string = sprintf('"%s" -netlist "%s%s_opt.asc"',spicePath, filePath, fileName);
fprintf('Issuing command to write new netlist from optimized schematic\n%s\n',string);
system(string);
pause(0.1);


% run sim on _Opt schematic
fprintf('Issuing command to run LTspice sim on optmized netlist\n%s\n',RunLTstring_opt);
[status,cmdout] = system(RunLTstring_opt);
pause(LTSpice_simTime);
result = LTspice2Matlab(LTSpice_outputfile_opt); % parse output file
pause(0.1)

% get opt response, with quantization
for i = 1:result.num_variables
    if strcmp(result.variable_name_list{i},LTSpice_output_node)
        fresp_opt_quant = abs(result.variable_mat(i,:));
        phase_opt_quant = unwrap(angle(result.variable_mat(i,:)));

    end
end


% plot the target, optimized, and quantized optimized responses
if matchMode==1
    figure;
    semilogx(freqx,20*log10(target),'g');
    hold on;
    semilogx(freqx,20*log10(fresp_opt),'r');
    hold on;
    semilogx(freqx,20*log10(fresp_opt_quant),'b');
    title('Ampl Response Opt results');
    legend('target (dB)','opt (dB)','opt quant (dB)');
end

if matchMode==2 % phase only
    figure;
    semilogx(freqx,target,'g');
    hold on;
    semilogx(freqx,phase_opt,'r');
    hold on;
    semilogx(freqx,phase_opt_quant,'b');
    title('Phase Opt results (radians)');
    legend('target','opt','opt quant');
end

if matchMode==3 % phase and ampl
    figure;
    semilogx(freqx,20*log10(target(1:end/2)),'g');
    hold on;
    semilogx(freqx,20*log10(fresp_opt),'r');
    hold on;
    semilogx(freqx,20*log10(fresp_opt_quant),'b');
    title('Ampl Response Opt results');
    legend('target (dB)','opt (dB)','opt quant (dB)');

    figure;
    semilogx(freqx,target(end/2+1:end),'g');
    hold on;
    semilogx(freqx,phase_opt,'r');
    hold on;
    semilogx(freqx,phase_opt_quant,'b');
    title('Phase Response Opt results (radians)');
    legend('target','opt','opt quant');
end

fprintf('\n*******\n***** DONE! ******\n*******\n');
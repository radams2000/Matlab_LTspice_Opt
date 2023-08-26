function err = optLTspice(optParams)
% This is the function called by lsqnonlin.
% It modifies the stored-in-memory netlist to update the variable components
% (identified in "simControl.m"), and re-writes the netlist file.
% Then it kicks off a simulation, gets the results, and computes the error
% between the target and the computed frequency response.
% It may take several hundred passes (or more) before the optimizer decides
% it's done the best it can do.
% There are some try/catch loops in order to recover from occasional failed
% simulations that can be caused by external things that cause slow
% execution, such as cloud-based automatic backup or some unknown mystery
% of the OS.

global passCell;
numOptd = passCell{6};
OptLine=passCell{7} ;
nomParams=passCell{8};
netlist=passCell{9};
netlist_fname=passCell{10};
numlines_netlist=passCell{11};
RunLTstring=passCell{12};
LTSpice_outputfile=passCell{13};
target=passCell{5};
errWeights = passCell{3};
LTSpice_simTime = passCell{16};
LTSpice_output_node = passCell{17};
matchMode = passCell{18};
filePath=passCell{2};
persistent simCount;
persistent simRestartCount;

if(isempty(simCount))
    simCount=0;
end

if(isempty(simRestartCount))
    simRestartCount=0;
end



for k=1:numOptd
    netlist{OptLine(k)}(4) = cellstr(sprintf('%2.12e ',optParams(k)*nomParams(k)));
end


fprintf('\ncurrent component values\n');
for k=1:numOptd
     fprintf('%s %2.12e\n',char(netlist{OptLine(k)}(1)),optParams(k)*nomParams(k));
   
end

fid_wr_netlist = fopen(netlist_fname,'W');

for k=1:numlines_netlist
    N=size(netlist{k},2);
    thisLine = netlist{k};
    for n=1:N
        fprintf(fid_wr_netlist,'%s ',thisLine{n});
    end
    fprintf(fid_wr_netlist,'\n ');
end
fclose(fid_wr_netlist);
pause(0.1);


% run the simulation

[status,simresult] = system(RunLTstring);
pause(LTSpice_simTime);

if(status) 
    fprintf('ERROR, LTSpice sim failed to run. Check setup or increase LTSpice_simTime\n ');
    return;
end

try
    result = LTspice2Matlab(LTSpice_outputfile); % parse output file
catch
    fprintf(' ltspice2matlab error, re-trying ...\n');
    fclose('all');
    simRestartCount = simRestartCount+1;
    [status,simresult] = system(RunLTstring);
    pause(5); % wait long enough for any file lock-ups to resolve (like dropbocx etc)
    result = LTspice2Matlab(LTSpice_outputfile); % parse output file
end

pause(0.1);



% get latest response
for i = 1:result.num_variables
    if strcmp(result.variable_name_list{i},LTSpice_output_node)
        fresp = abs(result.variable_mat(i,:));
        if(matchMode==2 || matchMode==3)
            phase = unwrap(angle(result.variable_mat(i,:)));
        end

    end
end




if matchMode==1 % ampl only
    optCurrent = fresp;
end
if matchMode==2 % phase only
    optCurrent = phase;
end
if matchMode==3 % ampl and phase; concatentae vectors, target must also be concatenated ampl and freq
    optCurrent = [fresp phase];
end


if(length(target) ~= length(optCurrent))
    fprintf('ERROR, something went wrong with the LTSpice sim. Re-trying...\n');
    fclose('all');
    simRestartCount = simRestartCount+1;
    [status,simresult] = system(RunLTstring);
    pause(5); % wait long enough for any file lock-ups to resolve (like dropbocx etc)
    result = LTspice2Matlab(LTSpice_outputfile); % parse output file
    % get opt response
    for i = 1:result.num_variables
        if strcmp(result.variable_name_list{i},LTSpice_output_node)
            fresp = abs(result.variable_mat(i,:));
            if(matchMode==2 || matchMode==3)
                phase = unwrap(angle(result.variable_mat(i,:)));
            end
        end
    end
end


err = target - optCurrent; % this is the error between target and the current freq response

err = err.*errWeights; % for frequency-dependent optimization, set in setTarget.m

fprintf('\ncurrent rms freq resp error = %2.6e\n',rms(err));
fclose('all');
simCount = simCount+1;

fprintf('sim count, sim restart count = %d %d\n',simCount,simRestartCount);


end
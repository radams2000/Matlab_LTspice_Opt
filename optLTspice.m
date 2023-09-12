function err = optLTspice(optParams)
% This is the function called by lsqnonlin.
% It modifies the stored-in-memory netlist to update the variable components
% (identified in "simControl.m"), and re-writes the netlist file.
% Then it kicks off a simulation, gets the results, and computes the error
% between the target and the computed frequency response.
% It may take several hundred passes (or more) before the optimizer decides
% it's done the best it can do.


global passCell;

persistent simCount;


if(isempty(simCount))
    simCount=0;
end


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
LTSpice_output_node = passCell{17};
matchMode = passCell{18};





for k=1:numOptd
    % netlist{OptLine(k)}(4) = cellstr(sprintf('%2.12e ',optParams(k)*nomParams(k)));
     netlist{OptLine(k)}(4) = cellstr(sprintf('%2.12e ',nomParams(k)*exp(optParams(k)) ));
end


fprintf('\ncurrent component values\n');
for k=1:numOptd
    fprintf('%s %2.12e\n',char(netlist{OptLine(k)}(1)),nomParams(k)*exp(optParams(k)));

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
pause(0.05);


% run the simulation

[status,simresult] = system(RunLTstring);

pause(0.1);

if(status)
    fprintf('ERROR, LTSpice sim failed to run');
    return;
end

result = LTspice2Matlab(LTSpice_outputfile); % parse output file

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

err = target - optCurrent; % this is the error between target and the current freq response

err = err.*errWeights; % for frequency-dependent optimization, set in setTarget.m

fprintf('\ncurrent rms freq resp error = %2.6e\n',rms(err));
fclose('all');
simCount = simCount+1;

fprintf('sim count = %d\n',simCount);


end
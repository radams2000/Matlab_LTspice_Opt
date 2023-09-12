function [] = update_schematic(passCell,simControl)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%global passCell;
numOptd = passCell{6};
OptLine=passCell{7} ;
nomParams=passCell{8};
netlist=passCell{9};
filePath=passCell{2};
fileName=passCell{4};
X=passCell{14};
simControlInstTol=simControl{7};
simControlOptInstNames=simControl{4};



% read in schematic to update
schem_fname = sprintf('%s%s.asc',filePath,fileName); % schem filename
fid = fopen(schem_fname);
schem = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
schem = schem{1};
schem = regexp(schem,' ','split'); % schem(n) contains the entire line; and each line contains a cell array, with every entry in a different cell
numlines_schem = size(schem,1);

changeNext=0;
roundStringNext = 'E96';
for k=1:numlines_schem
    if(changeNext==1) 
        newVal = round63(instValNext,roundStringNext);
        schem{k}(3) = cellstr(sprintf('%1.3e',newVal));
        fprintf('Inst, opt val, quantized val = %s %1.3e %1.3e\n',char(instNm),instValNext,newVal);
    end

    changeNext=0;

    if(strcmp(char(schem{k}(2)),'InstName'))

        instNm=(schem{k}(3));
        % scan all the optimized instances to see if this instance needs its value updated
        changeNext=0;
        for kk=1:numOptd    
            if(strcmp((netlist{OptLine(kk)}(1)),instNm))
                % find the index to this instance in
                % simControlOPtInstNames so that we know which tolerance to use
                for j=1:length(simControlOptInstNames)
                    if strcmp(char(instNm),char(simControlOptInstNames(j)))
                        xx = j;
                    end
                end
                % next line has the value to change
                changeNext=1;
                % instValNext = X(kk)*nomParams(kk);
                instValNext = nomParams(kk)*exp(X(kk));
                roundStringNext = char(simControlInstTol(xx));
            end
        end
    end
%fprintf('%s\n',char(line));
end

% write new schem file
schem_fname = sprintf('%s%s_opt.asc',filePath,fileName); % netlist filename

fid_wr_schem = fopen(schem_fname,'W');
for k=1:numlines_schem
    N=size(schem{k},2);
    thisLine = schem{k};
    for n=1:N
        fprintf(fid_wr_schem,'%s ',thisLine{n});
    end
    fprintf(fid_wr_schem,'\n');
end
fclose(fid_wr_schem);

end
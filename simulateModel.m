function result = simulateModel(spicePath, fileName, filePath)
% SIMULATEMODEL is used to simulate a spice netlist. It takes the following
% Inputs: 
% spicePath: The path to the LT SPice installation start file
% for example 'C:\Program Files\LTC\LTspiceXVII\start.exe' This should be a
% string
% filename: Name of the netlist file to be simulated. Do not use .net
% extension. This can be a string / char or even an int or float
% filepath: Provide the complete file path to the netlist and end it with a
% backslash (\). This has to be a string.
if ischar(fileName)
    filename = fileName;
else
    filename = num2str(fileName);
end


% write netlist
string = sprintf('"%s" -netlist "%s%s.asc"',spicePath, filePath, filename);
fprintf('%s\n',string);

dos(string);



string = sprintf('start "LTspice" "%s" -b "%s%s.net"',spicePath, filePath, filename);
fprintf('%s\n',string);

dos(string);

outputfile = sprintf('%s%s.raw', filePath,filename);
fprintf('output filename = %s\n',outputfile);

pause(10); %This 5 seconds pause is needed so that the the LT spice analysis is completed and results are available to read. 
          %If it is a very big analysis you may need to increse from 5s to higher time.

result = LTspice2Matlab(outputfile);

end
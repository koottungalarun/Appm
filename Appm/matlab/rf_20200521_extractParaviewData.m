% Example for extracting data from Paraview Output file of a point over
% time

filename = "Testcase-Euler-fixed-Efield0.csv";
data = importdata(filename);
idx(1) = find(strcmp(data.colheaders, """avg(Ions velocity (2))"""));
idx(2) = find(strcmp(data.colheaders, """avg(Electrons velocity (2))"""));
idx(3) = find(strcmp(data.colheaders, """avg(Ions pressure)"""));
idx(4) = find(strcmp(data.colheaders, """avg(Electrons pressure)"""));

clear time
% Read time values from HDF5 files, or from text file
isReadTimeFromTextFile = true;
if isReadTimeFromTextFile
    temp = importdata('timesteps.dat');
    time = temp.data(:,2);
else
%     % Read from HDF5 files:
%     time = [];
%     for i = 0 : 10
%         filename = sprintf('appm-%d.h5', i);
%         time(end+1) = h5read(filename, '/time');
%     end
end


plot(time, data.data(:,idx))
grid on

legnd = {"Ions velocity"
    "Electrons velocity"
    "Ions pressure"
    "Electrons pressure"};

legend(legnd, 'Location', 'NW')

function data = loadPLAfiles()

[filename, pathname] = uigetfile( ...
{  '*.jpg;*.jpeg;*.tif;*.png;*.gif;*.bmp','All Image Files';}, ...
   'Pick a PLA files', ...
   'MultiSelect', 'on');

if isequal(filename,0)
    error('No files are loaded.');
end

if iscell(filename)==0 %cell check
    data = imread([pathname filename]); %load image
else
    n = length(filename);
    data = struct([]); %memory allocation
    for i = 1:n
        data{i} = imread([pathname filename{i}]); %load image
    end
end
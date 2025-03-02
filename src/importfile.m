function data = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   data = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   data = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   data = importfile('energydata_complete.csv', 1, 19736);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/03/05 14:25:15

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]);
rawStringColumns = string(raw(:, 1));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = table;
data.date = rawStringColumns(:, 1);
data.Appliances = cell2mat(rawNumericColumns(:, 1));
data.lights = cell2mat(rawNumericColumns(:, 2));
data.T1 = cell2mat(rawNumericColumns(:, 3));
data.RH_1 = cell2mat(rawNumericColumns(:, 4));
data.T2 = cell2mat(rawNumericColumns(:, 5));
data.RH_2 = cell2mat(rawNumericColumns(:, 6));
data.T3 = cell2mat(rawNumericColumns(:, 7));
data.RH_3 = cell2mat(rawNumericColumns(:, 8));
data.T4 = cell2mat(rawNumericColumns(:, 9));
data.RH_4 = cell2mat(rawNumericColumns(:, 10));
data.T5 = cell2mat(rawNumericColumns(:, 11));
data.RH_5 = cell2mat(rawNumericColumns(:, 12));
data.T6 = cell2mat(rawNumericColumns(:, 13));
data.RH_6 = cell2mat(rawNumericColumns(:, 14));
data.T7 = cell2mat(rawNumericColumns(:, 15));
data.RH_7 = cell2mat(rawNumericColumns(:, 16));
data.T8 = cell2mat(rawNumericColumns(:, 17));
data.RH_8 = cell2mat(rawNumericColumns(:, 18));
data.T9 = cell2mat(rawNumericColumns(:, 19));
data.RH_9 = cell2mat(rawNumericColumns(:, 20));
data.T_out = cell2mat(rawNumericColumns(:, 21));
data.Press_mm_hg = cell2mat(rawNumericColumns(:, 22));
data.RH_out = cell2mat(rawNumericColumns(:, 23));
data.Windspeed = cell2mat(rawNumericColumns(:, 24));
data.Visibility = cell2mat(rawNumericColumns(:, 25));
data.Tdewpoint = cell2mat(rawNumericColumns(:, 26));
data.rv1 = cell2mat(rawNumericColumns(:, 27));
data.rv2 = cell2mat(rawNumericColumns(:, 28));
data(1,:) = []; % delete first row

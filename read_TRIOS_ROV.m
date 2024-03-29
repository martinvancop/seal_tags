clear all; close all

indir = '/Users/ioulianikolskaia/Boulot/my_ENCADREMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/PS117_Compilation_Martin/PS117_ROV_Data'

%
% READ ROV TRIOS DATA
% (c) Martin V, 30 mars 2021
%

%--------------------------------------------------------------------------
% Scan list of directories
%--------------------------------------------------------------------------
file_list = ls(indir)
Filedir = strsplit(file_list)
N_files = size(Filedir,2) - 1

%--------------------------------------------------------------------------
% Read file
%--------------------------------------------------------------------------
delimiter = ',';

i_file = 4 % VA DE 1-10 ICI

    % Create full file path / file name
    indir_full = strcat( string(indir), '/', string(Filedir(i_file)) );
    zfiles = strsplit(ls( indir_full ));
    Filename(i_file) = zfiles(1);
    infile_full = strcat(indir_full,'/', string(Filename(i_file)));
    
    %---------------------
    % Read ancillary data
    %---------------------

    startRow = 3;

    % Read columns of data as text:
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(infile_full,'r');

    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error
    % occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    % Close the text file.
    fclose(fileID);

    % Convert the contents of columns containing numeric text to numbers.
    % Replace non-numeric text with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[2,3,5,6,7,8,9]
        % Converts text in the input cell array to numbers. Replaced non-numeric text with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1)
            % Create a regular expression to detect and remove non-numeric prefixes and suffixes.
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

    % Split data into numeric and string columns.
    rawNumericColumns = raw(:, [2,3,5,6,7,8,9]);
    rawStringColumns = string(raw(:, [1,4,10,11,12,13]));

    % Make sure any text containing <undefined> is properly converted to an <undefined> categorical
    idx = (rawStringColumns(:, 2) == "<undefined>");
    rawStringColumns(idx, 2) = "";
    
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

    % Allocate imported array to column variable names
    ztime(:) = rawStringColumns(:, 1);
    TRIOS_ROV.PingDepth(i_file,:) = cell2mat(rawNumericColumns(:, 1));
    TRIOS_ROV.PressureDepth(i_file,:) = cell2mat(rawNumericColumns(:, 2));
    RadianceSensorSerial1 = categorical(rawStringColumns(:, 2));
    TRIOS_ROV.RadianceInclinationdegrees(i_file,:) = cell2mat(rawNumericColumns(:, 3));
    IrradianceSensorSerial1 = cell2mat(rawNumericColumns(:, 4));
    TRIOS_ROV.IrradianceInclinationdegrees(i_file,:) = cell2mat(rawNumericColumns(:, 5));
    TRIOS_ROV.lat(i_file,:) = cell2mat(rawNumericColumns(:, 6));
    TRIOS_ROV.lon(i_file,:) = cell2mat(rawNumericColumns(:, 7));
    GAPSTimestamp1 = rawStringColumns(:, 3);
    AutoFocusPicturePath = rawStringColumns(:, 4);
    ManualFocusPicturePath = rawStringColumns(:, 5);
    CAMTimestamp1 = rawStringColumns(:, 6);

    % Clear temporary variables
    clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns idx;
    
    % treat date array
    zsize = size(AutoFocusPicturePath,1); N_rec = zsize;
    
    for i_rec = 1:N_rec
        
       zz1 = strsplit( ztime(i_rec) , '-' );
       zyear = str2double(zz1(1));
       zmon = str2double(zz1(2));
       
       zz2 = strsplit( zz1(3), 'T' );
       zday = str2double(zz2(1));
       
       zz3 = strsplit( zz2(2), ':' );
       zhour = str2double(zz3(1));
       zmin = str2double(zz3(2));
       zsec = str2double(zz3(3));
       
       datearray(i_file,i_rec) = datetime(zyear,zmon,zday,zhour,zmin,zsec);

    end


    %------------------------
    % Read spectrometer data
    %------------------------
    
    % --- Read wavelengths
    fileID = fopen(infile_full,'r');
    startRow = 2;
    endRow = 2;
    formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    TRIOS_ROV.lambda(i_file,:) = [dataArray{1:end-1}];
    
    fclose(fileID);
    clearvars fileID startRow endRow formatSpec dataArray ans;

    % --- Read Radiances
    fileID = fopen(infile_full,'r');
    startRow = 3;
    formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    radiance = [dataArray{1:end-1}];
    
    TRIOS_ROV.radiance(i_file,:,:) = radiance(:,:);
    
    fclose(fileID);
    clearvars fileID startRow formatSpec dataArray ans;

    % --- Read Irradiances
    fileID = fopen(infile_full,'r');
    startRow = 3;
    formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    irradiance = [dataArray{1:end-1}];
    
    TRIOS_ROV.irradiance(i_file,:,:) = irradiance(:,:);

    fclose(fileID);
    clearvars fileID startRow formatSpec dataArray ans;
    
    clearvars radiance irradiance
    
    clearvars delimiter

%--------------------------------------------------------------------------
% QC Plot
%--------------------------------------------------------------------------

% --- Station name ---

figure; set(gcf, 'Position',  [100, 100, 800, 800])

subplot(4,1,1), hold on; box on;
plot(datearray(i_file,:), TRIOS_ROV.PingDepth(i_file,:), 'LineWidth', 2);
plot(datearray(i_file,:), TRIOS_ROV.PressureDepth(i_file,:), 'LineWidth', 2);
plot(datearray(i_file,:), TRIOS_ROV.PressureDepth(i_file,:) - TRIOS_ROV.PingDepth(i_file,:), 'LineWidth', 2);
legend('Ping', 'Pressure', 'Ice Draft')
ylabel('Depth (m)') 
title(strcat("TRIOS - ", cellstr(Filedir(i_file))), 'Interpreter', 'none');
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(4,1,2), hold on; box on;
plot(datearray(i_file,:), TRIOS_ROV.RadianceInclinationdegrees(i_file,:), 'LineWidth', 2), 
plot(datearray(i_file,:), TRIOS_ROV.IrradianceInclinationdegrees(i_file,:), 'LineWidth', 2), 
ylabel('Incl (�)'),
legend('Rad', 'Irr')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(4,1,3), hold on; box on; 
ylabel('mW/(m^2 nm sr)')
plot(datearray(i_file,:), max(TRIOS_ROV.radiance(i_file,:,:), [], 3), 'LineWidth', 2);
plot(datearray(i_file,:), nanmean(TRIOS_ROV.radiance(i_file,:,:),3), 'LineWidth', 2),
legend('Rad max', 'Rad mean')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(4,1,4), hold on; box on; 
ylabel(' mW/(m^2 nm)')
plot(datearray(i_file,:), max(TRIOS_ROV.irradiance(i_file,:,:), [], 3), 'LineWidth', 2);
plot(datearray(i_file,:), nanmean(TRIOS_ROV.irradiance(i_file,:,:),3), 'LineWidth', 2),
legend('Irrad max', 'Irrad mean')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

figure; hold on; box on;
plot(TRIOS_ROV.lon(i_file,:), TRIOS_ROV.lat(i_file,:), 'ko' ), xlabel('lon (�)'), ylabel('lat (�)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);


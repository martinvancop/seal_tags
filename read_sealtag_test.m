clear all, close all

%---------------------------------------------------------------------------
%
% read_sealtag_test : import data from seal tags and convert into .mat file
%
% start: Martin Vancoppenolle, May 2019
%
%---------------------------------------------------------------------------

% main script parameters
user   = 'Martin' % 'Martin' or 'Bastien'
i_read = 0      % 1 to read data (takes time!), 0 if not
i_save = 0      % 1 to save data into outfile
i_load = 1      % 0 to load data from outfile (to test the file)
i_plot = 1      % 0 to make a quick test plot

outfile = 'test.mat'

%-----------------
% Set directories
%-----------------

if ( strcmp(user,'Martin') == 1 )
   indir  = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/Tag_data_csv_PS117/SAV_CSV/';
   outdir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/Tag_data_csv_PS117/'
end

if ( strcmp(user,'Bastien') == 1 )
   indir   = 'ze_directory_of_bastien_for_input_files';
   outdir  = 'ze_directory_of_bastien_for_output_files';
end

%---------------------------------------
% Generate list of files Filename(1:Ns)
%---------------------------------------

file_list = ls(indir);
Filename = strsplit(file_list); % Filename is an array that gives the name of the file for each station
Ns = size(Filename,2) - 1;

%-----------------
% Encode metadata
%-----------------

% --- 'IceSt4_Heli1'
station = 'IceSt4_Heli1'; zfile_station = string([ '17A0530_' station '.csv' ]);
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190112 ; % YYYYMMDD
lat(zaddr)          = -70.38528 ; lon(zaddr) = -8.752222 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.48  % core length (m)
hs(zaddr)           = 0.007 % snow depth (m)
Ichla(zaddr)        = nan   % integrated chlorophyll (to be computed)

% --- 'ROV6'
station = 'ROV6'; zfile_station = string([ '17A0530_' station '.csv' ]);
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190117 ; % YYYYMMDD
lat(zaddr)          = -71.02668 ; lon(zaddr) = -13.4930235 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

% @Bastien - c'est ici que tu pourrais encoder le reste des metadonnees (a partir de ta
% super table), pour chacune des stations...

%--------------------------------
% Import tag data from csv files
%--------------------------------

if ( i_read == 1 );

    % Some generic instructions about file format
    delimiter = ',';
    startRow = 5;
    formatSpec = '%q%q%q%*q%*q%q%q%q%[^\n\r]';
    
    % Reading loop

    disp ({' ** ... Start reading LOOP (takes time)...'})

    for i=1:Ns

        i
        Percent_Complete=(i/Ns)*100
        zfile = [indir,'/',char(Filename(i))] % full path of the file

        % Open the text file.
        fileID = fopen(zfile,'r');

        % Read columns of data according to the format.
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

        for col=[1,2,3,4,5,6]
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

        % Exclude rows with blank cells
        I = any(cellfun(@(x) isempty(x) || (ischar(x) && all(x==' ')),raw),2); % Find row with blank cells
        raw(I,:) = [];

        % Allocate imported array to column variable names
        zDate      = cell2mat(raw(:, 1));
        zDepth     = cell2mat(raw(:, 2));
        zT         = cell2mat(raw(:, 3));
        intaZ      = cell2mat(raw(:, 4));
        zLL1       = cell2mat(raw(:, 5));
        zLL2       = cell2mat(raw(:, 6));

        % Assing variables to exported arrays
        Nr(i,1)            = size(zDate,1);
        Date(i,1:Nr(i,1))  = zDate;
        Depth(i,1:Nr(i,1)) = zDepth;
        Temp(i,1:Nr(i,1))  = zT;
        LL1(i,1:Nr(i,1))   = zLL1;
        LL2(i,1:Nr(i,1))   = zLL2;

        % Clear temporary variables
        clearvars filename fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp I;
        clearvars zDate zDepth zT zLL1 zLL2 zfile zfile_station

    end

end

%------------------------
% Save into a .mat file
%------------------------

if ( i_save == 1 );
    
    save( [outdir, outfile], 'Ns', 'station', 'cruise_name', 'datestr', 'lat', 'lon', 'platform', 'cl', 'hs', 'Ichla', ...
        'Nr', 'Date', 'Depth', 'Temp', 'LL1', 'LL2' );
    
end

%------------------------
% Load .mat file
%------------------------

if ( i_load == 1 );
    
    test = load( [outdir, outfile] );
    
end

%-----------
% Test plot
%-----------

if ( i_plot == 1 ) 
    
    if ( ( i_save == 0 ) & ( i_load == 1 ) )
        LL1 = test.LL1;
        LL2 = test.LL2;
        Depth = test.Depth;
        Temp  = test.Temp;
    end
    
    figure; hold on;
    zaddr = find(~isnan(LL1) & LL1 ~= 0 & ~isnan(LL2) & LL2 ~= 0 );

    subplot(2,2,1), plot(Depth(zaddr)), ylabel('Depth');

    subplot(2,2,2), hold on; 
    plot(LL1(zaddr)), plot(LL2), ylabel('LL1,LL2')

    subplot(2,2,3), hold on;
    plot(Temp(zaddr)), ylabel('temperature')
    
end

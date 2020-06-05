clear all, close all

%---------------------------------------------------------------------------
%
% read_TAG : import data from seal tags and convert into .mat file
%
% start: Martin Vancoppenolle, May 2019
%
%---------------------------------------------------------------------------

% main script parameters
user   = 'Martin' % 'Martin' or 'Bastien'
i_read = 0     % 1 to read data (takes time!), 0 if not
i_save = 0     % 1 to save data into outfile
i_load = 1      % 0 to load data from outfile (to test the file)
i_plot = 1      % 0 to make a quick test plot

outfile = 'PSS117_TAG.mat'

%-----------------
% Set directories
%-----------------

if ( strcmp(user,'Martin') == 1 )
   indir  = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/Tag_data_csv_PS117/CORRECTED_DATA/';
   outdir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/MAT_FILES/'
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

%----------
% 'IceSt1'  (pas de donn√©es tag pour cette station)
%----------
% station = 'IceSt1'; zfile_station = string([ '17A0530_' station '.csv' ]);
% zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)
% 
% % add metadata
% cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
% datestr(zaddr)      = 20190103 ; % YYYYMMDD
% lat(zaddr)          = -69.06889 ; lon(zaddr) = -0.09944444 % latitude & longitude
% platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed
% 
% cl(zaddr)           = 1.18   % core length (m) sur un seul site
% hs(zaddr)           = 0.47  % snow depth (m) sur un seul site
% Ichla(zaddr)        = nan      % integrated chlorophyll (to be computed)  


%--------------
% 'IceSt2' (ce fichier contient les donnÈes du tag pour les stations Ice 2)
%--------------
station = 'IceSt2'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190110 ; % YYYYMMDD
lat(zaddr)          = -70.28556 ; lon(zaddr) = -9.985833 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed
% platform(zaddr)   = "ROV"

cl(zaddr)           = 0.9325   % core length moyenn√©e (m) sur 4 sites 
hs(zaddr)           = 0.03025  % snow depth moyenn√©e (m) sur 4 sites
Ichla(zaddr)        = nan      % integrated chlorophyll (to be computed)   
ecart_type(zaddr)   = 0.087    % STD on core length (m)

%-----------
% 'IceSt3'
%-----------
station = 'IceSt3'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190111 ; % YYYYMMDD
lat(zaddr)          = -70.33583 ; lon(zaddr) = -8.935278 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.67   % core length mean over 2 sites
hs(zaddr)           = 0.0378 % snow depth mean over 2 sites
Ichla(zaddr)        = nan    % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0.07   % STD on core length (m)

%--------------------
%   'IceSt4_Heli1'
%--------------------
station = 'IceSt4_Heli1'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190112 ; % YYYYMMDD
lat(zaddr)          = -70.38528 ; lon(zaddr) = -8.752222 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.72   % core length mean over 2 sites
hs(zaddr)           = 0.0185 % snow depth mean over 2 sites
Ichla(zaddr)        = nan    % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0.24   % STD on core length (m)

%----------------
% 'IceSt5_Heli2'
%----------------
station = 'IceSt5_Heli2'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190113 ; % YYYYMMDD
lat(zaddr)          = -70.50056 ; lon(zaddr) = -8.185 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.095   % core length mean over 2 sites
hs(zaddr)           = 0.0615  % snow depth mean over 2 sites
Ichla(zaddr)        = nan     % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0.125   % STD on core length (m)


%----------------
% 'IceSt6_Heli3'
%----------------
station = 'IceSt6_Heli3'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190116 ; % YYYYMMDD
lat(zaddr)          = -70.50111 ; lon(zaddr) = -9.166667 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.11   % core length mean over 2 sites
hs(zaddr)           = 0.186  % snow depth mean over 2 sites
Ichla(zaddr)        = nan    % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0.06   % STD on core length (m)


%-----------------
% 'IceSt7_Heli4'
%-----------------
station = 'IceSt7_Heli4'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190117 ; % YYYYMMDD
lat(zaddr)          = -70.63361 ; lon(zaddr) = -9.901667 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 0.82    % core length mean over 2 sites
hs(zaddr)           = 0.0905  % snow depth mean over 2 sites
Ichla(zaddr)        = nan     % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0.03    % STD on core length (m)


%-----------------
% 'IceSt8_Heli5'
%-----------------
station = 'IceSt8_Heli5'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190126 ; % YYYYMMDD
lat(zaddr)          = -65.20111 ; lon(zaddr) = -45.935 % latitude & longitude
platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = 1.88    % core length (m) sur un seul site
hs(zaddr)           = 0.122   % snow depth (m) sur un seul site
Ichla(zaddr)        = nan     % integrated chlorophyll (to be computed)
ecart_type(zaddr)   = 0       % Ècart type de la longueur des carottes de glace                    


%-----------
%  'ROV1'
%-----------
station = 'ROV1'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190103 ; % YYYYMMDD
lat(zaddr)          = -69.1024863 ; lon(zaddr) = -0.0379608 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)


%-----------
%  'ROV 2'
%-----------
station = 'ROV2'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190109 ; % YYYYMMDD
lat(zaddr)          = -70.1951665 ; lon(zaddr) = -11.248467 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)


%----------
% 'ROV 3'
%----------
station = 'ROV3'; zfile_station = string([ '17A0530_' station '.csv' ]);
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190110 ; % YYYYMMDD
lat(zaddr)          = -70.2991597 ; lon(zaddr) = -10.0443703 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)


%----------
% 'ROV 4'   (pas de donn√©es tag pour la station ROV 4)
%----------
% station = 'ROV4'; zfile_station = string([ '17A0530_' station '.csv' ]);
% zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)
% 
% % add metadata
% cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
% datestr(zaddr)      = 20190111 ; % YYYYMMDD
% lat(zaddr)          = -70.3516683 ; lon(zaddr) = -8.975854 % latitude & longitude
% platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed
% 
% cl(zaddr)           = nan;  % core length (m)
% hs(zaddr)           = nan;  % snow depth (m)
% Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)


%----------
% 'ROV 5'
%----------
station = 'ROV5'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190117 ; % YYYYMMDD
lat(zaddr)          = -70.4590423 ; lon(zaddr) = -8.3512282 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)


%-----------
%  'ROV6'
%-----------
station = 'ROV6'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190118 ; % YYYYMMDD
lat(zaddr)          = -71.02668 ; lon(zaddr) = -13.4930235 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

%---------
% 'ROV 7'
%---------
station = 'ROV7'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190119 ; % YYYYMMDD
lat(zaddr)          = -71.2415445 ; lon(zaddr) = -19.6462462 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

%----------
% 'ROV 8'
%----------
station = 'ROV8'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190125 ; % YYYYMMDD
lat(zaddr)          = -67.6691357 ; lon(zaddr) = -46.3136623 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

%----------
% 'ROV 9'
%----------
station = 'ROV9'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190125 ; % YYYYMMDD
lat(zaddr)          = -67.5325338 ; lon(zaddr) = -46.3817188 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

%----------
% 'ROV 10'
%----------
station = 'ROV10'; zfile_station = string([ '17A0530_' station '.csv' ])
zaddr   = find( strcmp(string(Filename),zfile_station ) == 1); station_name(zaddr) = string(station);

% add metadata
cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
datestr(zaddr)      = 20190126 ; % YYYYMMDD
lat(zaddr)          = -65.5224318 ; lon(zaddr) = -46.113281 % latitude & longitude
platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed

cl(zaddr)           = nan;  % core length (m)
hs(zaddr)           = nan;  % snow depth (m)
Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)

% %????????????
% 
% % --- 'IceSt4_Heli1'
% station = 'IceSt4_Heli1'; zfile_station = string([ '17A0530_' station '.csv' ]);
% zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)
% 
% % add metadata
% cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
% datestr(zaddr)      = 20190112 ; % YYYYMMDD
% lat(zaddr)          = -70.38528 ; lon(zaddr) = -8.752222; % latitude & longitude
% platform(zaddr)     = "LARM"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed
% 
% cl(zaddr)           = 1.48  % core length (m)
% hs(zaddr)           = 0.007 % snow depth (m)
% Ichla(zaddr)        = nan   % integrated chlorophyll (to be computed)
% 
% % --- 'ROV6'
% station = 'ROV6'; zfile_station = string([ '17A0530_' station '.csv' ]);
% zaddr   = find( strcmp(string(Filename),zfile_station ) == 1)
% 
% % add metadata
% cruise_name(zaddr)  = "PS117" ; % "PSS117" or "OPTI2018" or ...
% datestr(zaddr)      = 20190117 ; % YYYYMMDD
% lat(zaddr)          = -71.02668 ; lon(zaddr) = -13.4930235; % latitude & longitude
% platform(zaddr)     = "ROV"  ; % "LARM" or "ROV" or "ICET" depending on which platform the tag was deployed
% 
% cl(zaddr)           = nan;  % core length (m)
% hs(zaddr)           = nan;  % snow depth (m)
% Ichla(zaddr)        = nan;  % integrated chlorophyll (to be computed)
% 
% % @Bastien - c'est ici que tu pourrais encoder le reste des metadonnees (a partir de ta
% % super table), pour chacune des stations...

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
        DateNum(i,1:Nr(i,1))  = zDate - 1.; % there is an 1 day lag for tag dates
        DateNum(i,1:Nr(i,1))  = addtodate( DateNum(i,1:Nr(i,1)), 1900, 'y' );
        DateTime(i,1:Nr(i,1)) = datetime( DateNum(i,1:Nr(i,1)), 'ConvertFrom', 'datenum' );
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
    
    save( [outdir, outfile], 'Ns', 'cruise_name', 'station_name', 'datestr', 'lat', 'lon', 'platform', 'cl', 'hs', 'Ichla', ...
        'Nr', 'DateNum', 'DateTime', 'Depth', 'Temp', 'LL1', 'LL2' );
    
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
        LL1          = test.LL1;
        LL2          = test.LL2;
        Depth        = test.Depth;
        Temp         = test.Temp;
        DateNum      = test.DateNum;
        Nr           = test.Nr;
        station_name = test.station_name
    end
    
    figure; hold on; % donnÈes brutes
    
    zaddr = find(~isnan(LL1) & LL1 ~= 0 & ~isnan(LL2) & LL2 ~= 0 );

    subplot(2,2,1), plot(Depth(zaddr)), ylabel('Depth');

    subplot(2,2,2), hold on; 
    plot(LL1(zaddr)), plot(LL2), ylabel('LL1,LL2')

    subplot(2,2,3), hold on;
    plot(Temp(zaddr)), ylabel('temperature')
    
    subplot(2,2,4), hold on;
    plot(DateNum(zaddr),Depth(zaddr), 'k.'), datetick
    
    figure; hold on %--- histograms
    subplot(2,2,1); hold on; box on;
    histogram(Depth(zaddr)), ylabel('Depth(m)')
    
    subplot(2,2,2); hold on; box on;
    histogram(Temp(zaddr)), ylabel('T (∞C)')
    
    subplot(2,2,3); hold on; box on;
    histogram(Temp(zaddr)), ylabel('LL1')
    
    subplot(2,2,4); hold on; box on;
    histogram(Temp(zaddr)), ylabel('LL2')
    
    figure; hold on %--- time series
    i_sta = 1; % retained station for ze plot

    subplot(3,1,1); hold on; box on;
    plot( DateNum(i_sta,1:Nr(i_sta)), -Depth(i_sta,1:Nr(i_sta)), 'k.' ); datetick('x','dd/mm HH:MM'); 
    title(station_name(i_sta), 'Interpreter', 'none')
    ylabel("Depth (m)");
    
    subplot(3,1,2); hold on; box on;
    plot( DateNum(i_sta,1:Nr(i_sta)), Temp(i_sta,1:Nr(i_sta)), 'k.' ); datetick('x','dd/mm HH:MM'); 
    ylabel("Temperature (∞C)")
    
    subplot(3,1,3); hold on; box on;
    plot( DateNum(i_sta,1:Nr(i_sta)), LL1(i_sta,1:Nr(i_sta)), 'g.' ); datetick('x','dd/mm HH:MM'); 
    plot( DateNum(i_sta,1:Nr(i_sta)), LL2(i_sta,1:Nr(i_sta)), 'b.' ); datetick('x','dd/mm HH:MM'); 
    ylabel("LL1,LL2")
    
    figure; hold on %--- profiles
    subplot(1,2,1); hold on; box on;
    plot( Temp(i_sta,:), - Depth(i_sta,:), 'k.' ); title(station_name(i_sta), 'Interpreter', 'none')
    xlabel("Temperature (∞C)"); ylabel("Depth (m)")
    subplot(1,2,2); hold on; box on;
    plot( LL1(i_sta,:), - Depth(i_sta,:), 'g.' ); 
    plot( LL2(i_sta,:), - Depth(i_sta,:), 'b.' ); 
    xlabel("LL1,LL2"); ylabel("Depth (m)")
    
    figure; hold on %--- number of records per file
    plot(DateNum(:,1),Nr,'ksq', 'MarkerFaceColor', 'blue'); ylabel('Number of records per file'); datetick  
    
end

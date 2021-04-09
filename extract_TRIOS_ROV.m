clear all; close all

indir = '/Users/ioulianikolskaia/Boulot/my_ENCADREMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/PS117_Compilation_Martin/PS117_ROV_Data'

%
% READ ROV TRIOS DATA and store I1f, I2f, EPAR, QPAR
% (c) Martin V, 30 mars 2021
%

outfile = 'PSS117_TRIOS_ROV.mat'
outdir = '/Users/ioulianikolskaia/Boulot/my_ENCADREMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/MAT_FILES/'

i_save = 1

%--------------------------------------------------------------------------
% Scan list of directories
%--------------------------------------------------------------------------
file_list = ls(indir);
Filedir = strsplit(file_list)
N_files = size(Filedir,2) - 1

for i_file = 1:N_files

%==========================================================================
%
% Read data
%
%==========================================================================
    
    delimiter = ',';

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
    TRIOS_ROV.N_rec(i_file) = size(rawStringColumns,1)
    ztime(1:TRIOS_ROV.N_rec(i_file)) = rawStringColumns(:, 1);
    TRIOS_ROV.PingDepth(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 1));
    TRIOS_ROV.PressureDepth(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 2));
    RadianceSensorSerial1 = categorical(rawStringColumns(:, 2));
    TRIOS_ROV.RadianceInclinationdegrees(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 3));
    IrradianceSensorSerial1 = cell2mat(rawNumericColumns(:, 4));
    TRIOS_ROV.IrradianceInclinationdegrees(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 5));
    TRIOS_ROV.lat(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 6));
    TRIOS_ROV.lon(i_file,1:TRIOS_ROV.N_rec(i_file)) = cell2mat(rawNumericColumns(:, 7));
    GAPSTimestamp1 = rawStringColumns(:, 3);
    AutoFocusPicturePath = rawStringColumns(:, 4);
    ManualFocusPicturePath = rawStringColumns(:, 5);
    CAMTimestamp1 = rawStringColumns(:, 6);

    % Clear temporary variables
    clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns idx;
    
    % treat date array
    %zsize = size(AutoFocusPicturePath,1); %N_rec = zsize;
    
    for i_rec = 1:TRIOS_ROV.N_rec(i_file)
        
       zz1 = strsplit( ztime(i_rec) , '-' );
       zyear = str2double(zz1(1));
       zmon = str2double(zz1(2));
       
       zz2 = strsplit( zz1(3), 'T' );
       zday = str2double(zz2(1));
       
       zz3 = strsplit( zz2(2), ':' );
       zhour = str2double(zz3(1));
       zmin = str2double(zz3(2));
       zsec = str2double(zz3(3));
       
       TRIOS_ROV.date(i_file,i_rec) = datetime(zyear,zmon,zday,zhour,zmin,zsec);

    end
    
    clear ztime

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
    
    TRIOS_ROV.radiance(i_file,1:TRIOS_ROV.N_rec(i_file),:) = radiance(:,:);
    
    fclose(fileID);
    clearvars fileID startRow formatSpec dataArray ans;

    % --- Read Irradiances
    fileID = fopen(infile_full,'r');
    startRow = 3;
    formatSpec = '%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    irradiance = [dataArray{1:end-1}];
    
    TRIOS_ROV.irradiance(i_file,1:TRIOS_ROV.N_rec(i_file),:) = irradiance(:,:);
    TRIOS_ROV.irradiance( find( TRIOS_ROV.irradiance < 0. ) ) = 0.;

    fclose(fileID);
    clearvars fileID startRow formatSpec dataArray ans;
    
    clearvars radiance irradiance
    
    clearvars delimiter
    
    clear rawNumericColumns
    
    %--------------------------------------------------------------------------
    % Filter transmittances
    %--------------------------------------------------------------------------
    
    if ( i_file == 1 );

        load("filter_responses.mat",'OPTI_18_filter_raw');
        lambda_f = OPTI_18_filter_raw.lambda;
        t1f_raw  = OPTI_18_filter_raw.rsp_green;
        t2f_raw  = OPTI_18_filter_raw.rsp_blue;
        clear OPTI_18_filter_raw

        % --- Interpolate filter transmittances on TRIOS-wavelength, then normalise
        TRIOS_ROV.t1f  = interp1(lambda_f,t1f_raw,TRIOS_ROV.lambda,'spline', 0.);
        TRIOS_ROV.t2f  = interp1(lambda_f,t2f_raw,TRIOS_ROV.lambda,'spline', 0);

        TRIOS_ROV.t1f( find( TRIOS_ROV.t1f < 0. ) ) = 0.;
        TRIOS_ROV.t2f( find( TRIOS_ROV.t2f < 0. ) ) = 0.;

        dlambda = diff(TRIOS_ROV.lambda); dlambda(254)=dlambda(253); % wavelength differences in trios space

        TRIOS_ROV.t1f  = TRIOS_ROV.t1f ./ nansum( TRIOS_ROV.t1f .* dlambda ); % normalise transmittance functions to unity integral
        TRIOS_ROV.t2f  = TRIOS_ROV.t2f ./ nansum( TRIOS_ROV.t2f .* dlambda ); % normalise transmittance functions to unity integral

    end
  
%==========================================================================
%
% Processing
%
%==========================================================================

    %----------------------
    % Filtered Irradiances
    %----------------------
    % By filtered irradiances, I mean integrated irradiances assuming filter tranmsission

    % --- Load filter transmittances


    planck = 6.62e-34;
    clight = 3.e8;

    for i_rec = 1:TRIOS_ROV.N_rec(i_file)

        zIrr(1,:) = TRIOS_ROV.irradiance(i_file, i_rec, :);

        TRIOS_ROV.I1f(i_file,i_rec) = 1.0e-3 * nansum( zIrr .* TRIOS_ROV.t1f .* dlambda ); % W/m2/nm 
        TRIOS_ROV.I2f(i_file,i_rec) = 1.0e-3 * nansum( zIrr .* TRIOS_ROV.t2f .* dlambda ); 

        TRIOS_ROV.EPAR(i_file,i_rec) = 1.0e-3 * nansum( zIrr(29:117) .* dlambda(29:117)  ); % W/m2/nm

        TRIOS_ROV.QPAR(i_file,i_rec) = 1.0e-3 * nansum( zIrr(29:117) .* dlambda(29:117) .* TRIOS_ROV.lambda(29:117) ) ...
                           / planck / clight * 1.0e-9 / 6.02e18; % µE/m2/s

    end

end

% 
% % --- Linear regression coefficients for I1 vs I1f
% zI1 = reshape(I1,11*4*3,1); zI1f = reshape(I1f,11*4*3,1);
% zI2 = reshape(I2,11*4*3,1); zI2f = reshape(I2f,11*4*3,1);
% zaddr = find(~isnan(zI1) );
% ft1 = fittype({'x'});
% 
% [C1, G1] = fit(zI1(zaddr),zI1f(zaddr),ft1);
% [C2, G2] = fit(zI2(zaddr),zI2f(zaddr),ft1);

% imposer ordonnée à l'origine = 0

%==================================================================================================
% Save into a .mat file
%==================================================================================================

if ( i_save == 1 ) ;
    
   save( [outdir, outfile], 'TRIOS_ROV'  );
                     
end


%==========================================================================
%
% QC Plots
%
%==========================================================================

figure(1); set(gcf, 'Position',  [100, 100, 800, 300]); %I1f; I2f;

for i_file = 1:N_files
 
    subplot(2,5,i_file); hold on; box on;
    plot( TRIOS_ROV.date(i_file,:), TRIOS_ROV.I1f(i_file,:), 'g', 'LineWidth', 2 )
    plot( TRIOS_ROV.date(i_file,:), TRIOS_ROV.I2f(i_file,:), 'b', 'LineWidth', 2 )
    title(Filedir(i_file), 'Interpreter', 'none')
    ylim([0.0001, 0.5])
    set(gca, 'YScale', 'log')
    
    if (i_file == 1); legend('I^1_f', 'I^2_f' ); ylabel('I_f (W/m^2)'); end
    set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
    
end

figure(2); set(gcf, 'Position',  [100, 100, 800, 300]); %E_{PAR}%
for i_file = 1:N_files
 
    subplot(2,5,i_file); hold on; box on;
    plot( TRIOS_ROV.date(i_file,:), TRIOS_ROV.EPAR(i_file,:), 'b', 'LineWidth', 2 )
    title(Filedir(i_file), 'Interpreter', 'none')
    %ylim([0.0001, 0.5])
    %set(gca, 'YScale', 'log')
    
    if (i_file == 1); ylabel('E_{PAR} (W/m^2)'); end
    set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
    
end

figure(3); set(gcf, 'Position',  [100, 100, 800, 300]); %Q_{PAR}
for i_file = 1:N_files
 
    subplot(2,5,i_file); hold on; box on;
    plot( TRIOS_ROV.date(i_file,:), TRIOS_ROV.QPAR(i_file,:), 'b', 'LineWidth', 2 )
    title(Filedir(i_file), 'Interpreter', 'none')
    %ylim([0.0001, 0.5])
    %set(gca, 'YScale', 'log')
    
    if (i_file == 1); ylabel('Q_{PAR} (µE/m^2/s)'); end
    set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
    
end

figure(4); set(gcf, 'Position',  [100, 100, 500, 300]); %Q_{PAR}

subplot(1,3,1); box on;
plot( TRIOS_ROV.EPAR, TRIOS_ROV.QPAR, 'ksq'), ylabel('Q_{PAR} (µE/m^2/s)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
subplot(1,3,2); box on;
plot( TRIOS_ROV.EPAR, TRIOS_ROV.I1f, 'gsq'), ylabel('I^1_f (W/m^2)'), xlabel('E_{PAR} (W/m^2)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
subplot(1,3,3); box on; 
plot( TRIOS_ROV.EPAR, TRIOS_ROV.I2f, 'bsq'), ylabel('I^2_f (W/m^2)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
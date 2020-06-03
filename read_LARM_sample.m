clear all; close all;

%---------------------------------------------------------------------------
%
% read_LARM_sample : import data from L-ARM stations
%
% start: Martin Vancoppenolle, May 2020
%
%---------------------------------------------------------------------------
% STATION 2, HOLE 2, REPLICA B in the AIR is WRONG due to operator error

% PAR, transmittance, NDI ?
% add metadata (thickness, snow depth, ...)

% main script parameters

user = 'Martin'
i_read = 1;
i_plot = 1;
i_save = 1;
i_load = 0;

outfile = 'PSS117_TRIOS.mat'

%==================================================================================================
% Set directories, list of stations, list of files
%==================================================================================================

%-----------------
% Set directories
%-----------------

if ( strcmp(user,'Martin') == 1 )
   indir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/PS117_Compilation_Martin/PS117_Ice_Station_Data/PS117_L_Arm_Spectra_Klaus/'
   outdir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/Tag_data_csv_PS117/'
end
if ( strcmp(user,'Bastien') == 1 )
   indir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/PS117_Compilation_Martin/PS117_Ice_Station_Data/PS117_L_Arm_Spectra_Klaus/'
end

%---------------------------------------
% Hard-code list of stations
%---------------------------------------
N_sta = 11

station_name(1) = "St1"
station_name(2) = "St2"
station_name(3) = "St3"
station_name(4) = "St4"
station_name(5) = "St5"
station_name(6) = "St6"
station_name(7) = "St7"
station_name(8) = "St8"
station_name(9) = "Calib1"
station_name(10) = "Calib2"
station_name(11) = "Calib3"

%---------------------------------------
% Generate list of files Filename(1:Ns)
%---------------------------------------

file_list = ls(indir);
Filename = strsplit(file_list); % Filename is an array that gives the name of the file for each station
Nfiles = size(Filename,2) - 1;

%==================================================================================================
% Extract data from text files
%==================================================================================================

if ( i_read == 1 );

    %--- Generic information for file format 
    % (irradiance data)
    delimiter = '\t';
    startRow = 46;
    endRow = 300;
    formatSpec = '%f%f%[^\n\r]';
    
    % (date data)
    startRow_date   = 4;
    endRow_date     = 4;
    formatSpec_date = '%*s%{yyyy-MM-dd HH:mm:ss}D%[^\n\r]';
    
    %--- Create empty arrays
    N_hole(1:N_sta) = 0;
    Irr_air   = nan(N_sta,4,3,255);
    Irr_water = nan(N_sta,4,3,255);

    disp ({' ** ... Start reading LOOP (takes time)...'})
    for i=1:Nfiles %--- Main Loop (over files)
        
        i
        Percent_Complete=(i/Nfiles)*100
               
        %-----------------------------
        %--- Extract station metadata
        %-----------------------------
        % (find which sample the file corresponds to)

        zstation(i) = [ "aaa" ];
        zzz         = split( Filename(i), '_' );
        znstr       = size(zzz,1);
        
        if ( znstr == 3 ) %--- Calibration stations
           zcruise(i) = zzz(1);
           zsample(i) = zzz(2);
           zstation(i) = zzz(3);
           i_hole     = 1
           i_replica  = 1
        end
        if ( znstr == 4 ) %--- All other stations
            
           zcruise(i) = zzz(1);
           zsample(i) = zzz(4);
           zstation(i) = zzz(2)
           zhole(i)    = string(zzz(3))
           
           if ( strcmp(zhole(i), 'H1') == 1 ); i_hole = 1; end
           if ( strcmp(zhole(i), 'H2') == 1 ); i_hole = 2; end
           if ( strcmp(zhole(i), 'H3') == 1 ); i_hole = 3; end
           if ( strcmp(zhole(i), 'H4') == 1 ); i_hole = 4; end
           
        end
        
        % ( this is wrong to use i_sta as an index for array!!! )
        i_sta = find( strcmp(zstation(i), station_name ) == 1 ); % Station ID
        N_hole(i_sta) = max( [ N_hole(i_sta) i_hole ] );          % Hole ID

        if ( strcmp(zsample(i),'Air' ) == 1 ); c_medium(i) = "Air"; i_medium = 1; i_replica = 1; end
        if ( strcmp(zsample(i),'AirA' ) == 1 ); c_medium(i) = "Air"; i_medium = 1; i_replica = 1; end
        if ( strcmp(zsample(i),'AirB' ) == 1 ); c_medium(i) = "Air"; i_medium = 1; i_replica = 2; end
        if ( strcmp(zsample(i),'AirC' ) == 1 ); c_medium(i) = "Air"; i_medium = 1; i_replica = 3; end
        if ( strcmp(zsample(i),'WaterA' ) == 1 ); c_medium(i) = "Water"; i_medium = 2; i_replica = 1; end
        if ( strcmp(zsample(i),'WaterB' ) == 1 ); c_medium(i) = "Water"; i_medium = 2; i_replica = 2; end
        if ( strcmp(zsample(i),'WaterC' ) == 1 ); c_medium(i) = "Water"; i_medium = 2; i_replica = 3; end
        
        %----------------------
        % Read irradiance data
        %----------------------
        zfile = [indir,'/',char(Filename(i))]; % full path of the file

        fileID = fopen(zfile,'r');
        dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);

        zlambda = dataArray{:, 1};
        zIrr    = dataArray{:, 2};
        
        if ( i == 1 );
            lambda = zlambda;
        end
        if ( i_medium == 1 ) ;  Irr_air(i_sta,i_hole,i_replica,:) = zIrr; end
        if ( i_medium == 2 ) ;  Irr_water(i_sta,i_hole,i_replica,:) = zIrr; end

        clearvars zlambda zIrr
        
        %--------------------------------
        % Read sampling date & time
        %--------------------------------
        fileID = fopen(zfile,'r');
        dataArray = textscan(fileID, formatSpec_date, endRow_date-startRow_date+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow_date-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);

        date_sample(i_sta,i_hole,i_replica) = dataArray{:, 1};
        
        clearvars dataArray ans
        
    end
    
    % --- Clean irradiance arrays (remove negative values)
    Irr_air( find (Irr_air < 0.) ) = 0.; % remove negative values
    Irr_water( find (Irr_water < 0.) ) = 0.; % remove negative values
    
    clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
    
end

%==================================================================================================
% Extract diagnostics (delta & filtered irradiances)
%==================================================================================================

i2 = 41; % index of peak-wavelength of blue filter
i1 = 74; % index of peak-wavelength of green filter

%-------------------
% Delta Irradiances
%-------------------
% By delta-irradiance, I mean single wavelength irradiance, 
% taken at peak-wavelength of filter transmittance (W/m2/nm)
for i_sta = 1:N_sta; for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
   I1(i_sta,i_hole,i_rep) = Irr_water(i_sta,i_hole,i_rep,i1);
   I2(i_sta,i_hole,i_rep) = Irr_water(i_sta,i_hole,i_rep,i2);
end; end; end;

%----------------------
% Filtered Irradiances
%----------------------
% By filtered irradiances, I mean integrated irradiances assuming filter tranmsission

% --- Load filter transmittances
load("filter_responses.mat",'OPTI_18_filter_raw');
lambda_f = OPTI_18_filter_raw.lambda;
t1f_raw  = OPTI_18_filter_raw.rsp_green;
t2f_raw  = OPTI_18_filter_raw.rsp_blue;
clear OPTI_18_filter_raw

% --- Interpolate filter transmittances on TRIOS-wavelength, then normalise
t1f_trios  = interp1(lambda_f,t1f_raw,lambda,'spline', 0.);
t2f_trios  = interp1(lambda_f,t2f_raw,lambda,'spline', 0);

dlambda = diff(lambda); dlambda(255)=dlambda(254); % wavelength differences in trios space

t1f_trios  = t1f_trios ./ nansum(t1f_trios .* dlambda(:) ); % normalise transmittance functions to unity integral
t2f_trios  = t2f_trios ./ nansum(t2f_trios .* dlambda(:) ); % normalise transmittance functions to unity integral

% --- Integrate irradiances over all spectrum
for i_sta = 1:N_sta; for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
    zIrr(:) = Irr_water(i_sta,i_hole,i_rep,:);
    I1f(i_sta,i_hole,i_rep) = nansum( zIrr(:) .* t1f_trios(:) .* dlambda(:) );  
    I2f(i_sta,i_hole,i_rep) = nansum( zIrr(:) .* t2f_trios(:) .* dlambda(:) ); 
end; end; end;

% --- Linear regression coefficients for I1 vs I1f
zI1 = reshape(I1,11*4*3,1); zI1f = reshape(I1f,11*4*3,1);
zI2 = reshape(I2,11*4*3,1); zI2f = reshape(I2f,11*4*3,1);
zaddr = find(~isnan(zI1) );
ft1 = fittype({'x'});

[C1, G1] = fit(zI1(zaddr),zI1f(zaddr),ft1);
[C2, G2] = fit(zI2(zaddr),zI2f(zaddr),ft1);

% imposer ordonnée à l'origine = 0

%==================================================================================================
% Save into a .mat file
%==================================================================================================

if ( i_save == 1 );
    
save( [outdir, outfile], 'N_sta', 'N_hole', 'station_name', 'date_sample', 'I1', 'I2', 'I1f', 'I2f', 'lambda', 'Irr_air', 'Irr_water', ...
                         't1f_trios', 't2f_trios'  );
                     
    
end


%-----------
% Test plot
%-----------

if ( i_plot == 1 ) 
    
    % --- All irradiance data ---
    figure;
    subplot(2,1,1); hold on; box on;
    for i_sta = 1:N_sta; for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
       zIrr(:) = Irr_air(i_sta,i_hole,i_rep,:);
       plot(lambda, zIrr(:), 'LineWidth', 2)
    end; end; end;
            
    ylabel('Spectral Irradiance [mW/(m^2 nm)]')
    %xlabel('\lambda [nm]')
    title('All Ice station data (L-ARM)')
    
    set(gca, 'FontName', 'Helvetica LT Std')
    set(gca,'FontSize',12);
    set(gca, 'YScale', 'log')
    
    subplot(2,1,2); hold on; box on;
    for i_sta = 1:N_sta; for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
       zIrr(:) = Irr_water(i_sta,i_hole,i_rep,:);
       plot(lambda, zIrr(:), 'LineWidth', 2)
    end; end; end;
    %ylabel('Spectral Irradiance [mW/(m^2 nm)]')
    xlabel('\lambda [nm]')
    %title('All Ice station data (L-ARM)')
    
    set(gca, 'FontName', 'Helvetica LT Std')
    set(gca,'FontSize',12);
    set(gca, 'YScale', 'log')
    
    % irradiance per station ? AIR
    figure;
    for i_sta = 1:11
        subplot(3,4,i_sta); hold on; box on
        for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
            zIrr(:) = Irr_air(i_sta,i_hole,i_rep,:);
            plot(lambda, zIrr(:), 'LineWidth', 2)
        end; end;
        ylabel('Spectral Irradiance [mW/(m^2 nm)]')
        xlabel('\lambda [nm]')
        title(station_name(i_sta))
        ylim([1e-2,1e4])
    
        set(gca, 'FontName', 'Helvetica LT Std')
        set(gca,'FontSize',8);
        set(gca, 'YScale', 'log')
    end

    % irradiance per station - WATER
    figure;
    for i_sta = 1:11
        subplot(3,4,i_sta); hold on; box on
        for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
            zIrr(:) = Irr_water(i_sta,i_hole,i_rep,:);
            plot(lambda, zIrr(:), 'LineWidth', 2)
        end; end;
        ylabel('Spectral Irradiance [mW/(m^2 nm)]')
        xlabel('\lambda [nm]')
        title(station_name(i_sta))
        ylim([1e-5,1e2])
    
        set(gca, 'FontName', 'Helvetica LT Std')
        set(gca,'FontSize',8);
        set(gca, 'YScale', 'log')
    end
    
%     %----------------------------------------
%     figure % plot all spectra for station 2
%     %----------------------------------------
%     i = 0;
%     i_sta = 2
%     for i_hole = 1:4; for i_rep = 1:3
%         i = i+1;
%         subplot(4,3,i)
%         zIrr(:) = Irr_air(i_sta,i_hole,i_rep,:);
%         plot(lambda,zIrr)
%         title(strcat('Hole ', string(i_hole),  '; Rep ', string(i_rep)))
%         end
%     end
    
    %-------------------------------
    figure % Filter characteristics
    %-------------------------------
    subplot(2,1,1); hold on; 
    %plot(lambda_f,t1f_raw, 'g--', 'LineWidth', 1); plot(lambda_f,t2f_raw, 'b--', 'LineWidth', 1);
    p1 = plot(lambda,t1f_trios, 'g', 'LineWidth', 3); p2 =  plot(lambda,t2f_trios, 'b', 'LineWidth', 2);
    xlabel('\lambda'), ylabel('T_f'), title('Filter transmittance')
    legend( [p1, p2], { 'T^f_1', 'T^f_2' } )
    set(gca, 'FontName', 'Helvetica LT Std')
    set(gca,'FontSize',12);

    subplot(2,1,2); hold on; box on;
    for i_sta = 1:N_sta; for i_hole = 1:N_hole(i_sta); for i_rep = 1:3
       plot(I1(i_sta,i_hole,i_rep),I1f(i_sta,i_hole,i_rep),'gsq', 'MarkerFaceColor', 'g')
       plot(zI1(zaddr), C1.a * zI1(zaddr), 'g:')
       plot(I2(i_sta,i_hole,i_rep),I2f(i_sta,i_hole,i_rep),'bsq', 'MarkerFaceColor', 'b')
       plot(zI2(zaddr), C2.a * zI2(zaddr), 'b:')
    end; end; end;
    xlabel('Delta-irradiances (W/m^2/nm)'), ylabel('Filtered irradiances (W/m^2)'), title('Filter vs delta irradiances')
    
    set(gca, 'FontName', 'Helvetica LT Std')
    set(gca,'FontSize',12);

return
    
end
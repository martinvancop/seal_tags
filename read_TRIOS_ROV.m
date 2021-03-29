clear all; close all

indir = '/Users/ioulianikolskaia/Boulot/my_ENCADREMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/PS117/PS117_Compilation_Martin/Crowsnest'

%--------------------------------------------------------------------------
% Scan list of files
%--------------------------------------------------------------------------

file_list = ls(indir);
Filename = strsplit(file_list)
N_files = size(Filename,2) - 1

%--------------------------------------------------------------------------
% Read file
%--------------------------------------------------------------------------

i_file = 1 % VA DE 1-10 ICI

infile = strcat(string(indir),'/',string(Filename(i_file))); % generate file name
fID = fopen(infile);                                        % open file
[status,cmdout] = system(strcat({'wc -l '}, infile));       % count number of lines in file
N_lines(i_file) = cell2mat(textscan(cmdout,'%f'))           % store numbers of line in file

i_rec = 0;
i_read_meta = 0;
i_read_data = 0;

for i_lines = 1:N_lines(i_file)
   tline = fgetl(fID);
   
   % Detect beginning of spectrum data
   if ( isequal(tline,'[Spectrum]') ); 
       %disp(tline); 
       i_rec = i_rec + 1;  
       i_read_meta = 1;
       
   end
   
   % Extract data
   if ( i_read_meta == 1 )
       
       % Extract metadata and ancillary stuff
       ztext = split(tline);
       
       if ( isequal(string(ztext(1)),'DateTime' )); % Date
          zdate = string(ztext(3)); 
          zhour = string(ztext(4)); 
          
          zzz   = strsplit(zdate,'-');
          
          zyear = str2double(zzz(1));
          zmon  = str2double(zzz(2));
          zday  = str2double(zzz(3));
          
          zzz   = strsplit(zhour,':');
          zhour = str2double(zzz(1));
          zmin  = str2double(zzz(2));
          zsec  = str2double(zzz(3));
          
          datearray(i_rec) = datetime(zyear,zmon,zday,zhour,zmin,zsec);
          
       end
       
       if ( isequal(string(ztext(1)),'Pressure' ) )        ; zpres(i_rec) = str2double(cell2mat(ztext(3))); end
       if ( isequal(string(ztext(1)),'InclV' ) )           ; zinclv(i_rec) = str2double(cell2mat(ztext(3))); end
       if ( isequal(string(ztext(1)),'InclX' ) )           ; zinclx(i_rec) = str2double(cell2mat(ztext(3))); end
       if ( isequal(string(ztext(1)),'InclY' ) )           ; zincly(i_rec) = str2double(cell2mat(ztext(3))); end
       if ( isequal(string(ztext(1)),'IntegrationTime' ) ) ; zinteg_time(i_rec) = str2double(cell2mat(ztext(3))); end
       
       if ( isequal(tline,'[DATA]') ); i_read_data = 1; i_band = 1; end
       if ( isequal(tline,'[END] of [DATA]') ); i_read_data = 0; end
       
       % Extract spectra
       if ( ~isequal(tline,'[DATA]') & i_read_data == 1)
           zlambda(i_rec,i_band) = str2double(cell2mat(ztext(2)));
           zcounts(i_rec,i_band) = str2double(cell2mat(ztext(3)));
           i_band = i_band + 1;
           %return
       end

   end
       
   % Detect end of spectrum data
   if ( isequal(tline,'[END] of Spectrum') ); 
      %disp(tline); 
 
       i_read_meta = 0;
   end
      
end

% convert pressure into depth
zdepth = - zpres / ( 1025 * 9.81 ) * 100000.;

N_rec(i_file) = i_rec;

fclose(fID);

%--------------------------------------------------------------------------
% Treat date and time
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% QC plots
%--------------------------------------------------------------------------

figure; set(gcf, 'Position',  [100, 100, 800, 400])

subplot(3,3,1), plot(datearray,zpres), ylabel('P (dBar)'), title(strcat("ROV - ", cellstr(Filename(i_file))));
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,2), plot(datearray,zdepth), ylabel('z (m)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,3), plot(datearray,zinclv), ylabel('Incl_v')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,4), plot(datearray,zinclx), ylabel('Incl_x')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,5), plot(datearray,zincly), ylabel('Incl_y')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,6), plot(datearray,zinteg_time, 'k.'), ylabel('T_{integ}')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,7), hold on; ylabel('Intensity (mean and max), mW/(m^2 nm)')
plot(datearray,max(zcounts, [], 2),'b.');
plot(datearray,nanmean(zcounts,2),'k.'),
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

subplot(3,3,8), plot(zpres,max(zcounts, [], 2), 'k.'), ylabel('I_{max}, mW/(m^2 nm)'), xlabel('Pres (dBar)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);



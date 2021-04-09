clear all;
close all;

% ADJUST SCRIPT TO DO EVERYTHING
% LL1, LL2, PAR(LL1,LL2) vs EPAR
% https://fr.mathworks.com/help/matlab/data_analysis/linear-regression.html
% NEED ALSO TO DO RESIDUALS
%---------------------------------------------------------------------------
%
% calibrate_TAG_ROV : test calibration of TAG vs ROV data
%
%---------------------------------------------------------------------------

user = 'Martin' % 'Martin'

if ( strcmp(user,'Martin') == 1 )
   indir = '/Users/ioulianikolskaia/Boulot/my_ENCADREMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/MAT_FILES/'
end

if ( strcmp(user,'Bastien') == 1 )
   indir = '/home/utl1/Bureau/'
end

% Load tag and trios data

TAG       = load( [ indir, 'PSS117_TAG.mat' ] );
load( [ indir, 'PSS117_TRIOS_ROV.mat' ] );

% ----------
% POUR LL1 
% ----------

zI1f        = reshape( TRIOS_ROV.I1f, 10*99, 1 ); % déplie la matrice en un vecteur 1D

zDate_ROV   = reshape( TRIOS_ROV.date, 10*99, 1 ); % deplie

zLL1        = TAG.LL1;

zDate_TAG   = TAG.DateNum;

zLL1_ROV = nan(10,99); % le vecteur qui contient les données du TAG dans l'espace des données du TRIOS

for i_ROV = 1:10 % boucle sur tous les échantillonnages ROV
    
    for i_REC = 1:TRIOS_ROV.N_rec(i_ROV)
        
    %if ( zI1f(i_ROV) > 0 ); % ne travailler que sur les données non-nulles
    
        for i_TAG = 1:16;
       
            zaddr_valid = find( zLL1(i_TAG,:) > 0 ); % addresses des points valides 
            zdiff       = abs( 	datenum(TRIOS_ROV.date(i_ROV,i_REC))- zDate_TAG(i_TAG,zaddr_valid) ); % différences de dates
            zaddr_min   = find( zdiff < 1e-5); % find difference less than a minute
       
            if ( size(zaddr_min,2) > 0 )
       
                datetime(zDate_TAG(i_TAG,zaddr_min),'ConvertFrom','datenum');
                TRIOS_ROV.date(i_ROV,i_REC);
       
                zaddr_tag_station(i_ROV,i_REC) = i_TAG;
                zaddr_tag_index(i_ROV,i_REC)   = min(zaddr_min);
                zLL1_ROV(i_ROV,i_REC) = zLL1( i_TAG, zaddr_tag_index(i_ROV,i_REC) ); 
                % on a trouvé la donnée du TAG correspondant Ã  l'échantillonnage ROV
       
            end
        end
    end
end

TRIOS_valid = find (TRIOS_ROV.I1f > 1e-5);
figure; hold on; box on
plot (zLL1_ROV(TRIOS_valid),(log10(TRIOS_ROV.I1f(TRIOS_valid))), 'ksq', 'MarkerFaceColor', 'b');
ylabel('log10(TRIOS-ROV.I1f)')
xlabel('TAG.LL1')
% plot (log10(TRIOS_ROV.I1f), zLL1_ROV, 'ksq', 'MarkerFaceColor', 'b');

zaddr = find(~isnan(zLL1_ROV(TRIOS_valid)))
x = zLL1_ROV(TRIOS_valid);
y = log10(TRIOS_ROV.I1f(TRIOS_valid));
[p,s] = polyfit (x(zaddr), y(zaddr), 1)
plot ([500:900], p(2)+p(1)*[500:900], 'k:', 'LineWidth', 2)
err = std (y(zaddr)./x(zaddr)) % erreur pente
err1 = std(y(zaddr)) % erreur ordonnée Ã  l'origine
r = corrcoef (x(zaddr), y(zaddr)) % coefficient de corrélation linéaire
R = r.*r % coefficient de détermination

s.normr

set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

% Recherches (ne pas décommenter): 
% -------------------------------------------------------------------------------------------------------------------
% figure; hold on;
% % if ( log10(TRIOS_ROV.I1f) > -5 )
% % if ( TRIOS_ROV.I1f (10,99) > 1e-5 )
% if ((log10(TRIOS_ROV.I1f))~=-8.64)
% plot( log10(TRIOS_ROV.I1f), zLL1_ROV,'ksq','MarkerFaceColor', 'b');
% xlabel('log10(TRIOS_ROV.I1f)')
% ylabel('TAG.LL1')
% end

% for i_ROV = 1:10
%     for i_REC = 1:TRIOS_ROV.N_rec(i_ROV)
%     if ( log10(TRIOS_ROV.I1f(i_ROV,i_REC)) > -5 )
%         plot( log10(TRIOS_ROV.I1f), zLL1_ROV,'ksq','MarkerFaceColor', 'b');
%     end
%     end
%     
% end

% for i_ROV = 1:10
%     for i_REC = 1:TRIOS_ROV.N_rec(i_ROV)
%         TRIOS_valid = find (TRIOS_ROV.I1f(i_ROV,i_REC) > 1e-5);
%     end
% end
% figure; hold on;
% plot ((log10(TRIOS_ROV.I1f(TRIOS_valid))), zLL1_ROV(TRIOS_valid), 'ksq', 'MarkerFaceColor', 'b');


% figure; hold on;
% plot( log10(TRIOS_valid), zLL1_ROV,'ksq','MarkerFaceColor', 'b');
%  TRIOS_valid = find(TRIOS_ROV.I1f (:,:) > 0.005)
%  TRIOS_valid
% -------------------------------------------------------------------------------------------------------------------


% ----------
% POUR LL2 
% ----------


% zLL2        = TAG.LL2;
% 
% zDate_TAG   = TAG.DateNum;
% 
% zLL2_ROV = nan(10,99); % le vecteur qui contient les données du TAG dans l'espace des données du TRIOS
% 
% for i_ROV = 1:10 % boucle sur tous les échantillonnages ROV
%     
%     for i_REC = 1:TRIOS_ROV.N_rec(i_ROV)
%         
%         for i_TAG = 1:16;
%        
%             zaddr_valid = find( zLL2(i_TAG,:) > 0 ); % addresses des points valides 
%             zdiff       = abs( 	datenum(TRIOS_ROV.date(i_ROV,i_REC))- zDate_TAG(i_TAG,zaddr_valid) ); % différences de dates
%             zaddr_min   = find( zdiff < 1e-5) % find difference less than a minute
%        
%             if ( size(zaddr_min,2) > 0 )
%        
%                 datetime(zDate_TAG(i_TAG,zaddr_min),'ConvertFrom','datenum')
%                 TRIOS_ROV.date(i_ROV,i_REC)
%        
%                 zaddr_tag_station(i_ROV,i_REC) = i_TAG;
%                 zaddr_tag_index(i_ROV,i_REC)   = min(zaddr_min);
%                 zLL2_ROV(i_ROV,i_REC) = zLL2( i_TAG, zaddr_tag_index(i_ROV,i_REC) ); 
%                 % on a trouvé la donnée du TAG correspondant Ã  l'échantillonnage ROV
%        
%             end
%         end
%     end
% end

% TRIOS_valid = find (TRIOS_ROV.I2f > 1e-5);
% figure; hold on;
% plot ((log10(TRIOS_ROV.I2f(TRIOS_valid))), zLL2_ROV(TRIOS_valid), 'ksq', 'MarkerFaceColor', 'b');
% xlabel('log10(TRIOS-ROV.I2f)')
% ylabel('TAG.LL2')

% %plot (log10(TRIOS_ROV.I2f), zLL2_ROV, 'ksq', 'MarkerFaceColor', 'b');

% zaddr = find(~isnan(zLL2_ROV(TRIOS_valid)))
% x = log10(TRIOS_ROV.I2f(TRIOS_valid));
% y = zLL2_ROV(TRIOS_valid);
% p = polyfit (x(zaddr), y(zaddr), 1)
% plot ([-5:0], 59.1584*[-5:0]+795.8309, 'k:', 'LineWidth', 2)
% err = std (y(zaddr)/x(zaddr)) % erreur pente
% err1 = std(y(zaddr)) % erreur ordonnée Ã  l'origine
% r = corrcoef (x(zaddr), y(zaddr)) % coefficient de corrélation linéaire
% R = r*r % coefficient de détermination

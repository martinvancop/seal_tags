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

zLL         = TAG.LL1 + TAG.LL2;

zDate_TAG   = TAG.DateNum;

zLL_ROV     = nan(10,99); % le vecteur qui contient les données du TAG dans l'espace des données du TRIOS

for i_ROV = 1:10 % boucle sur tous les échantillonnages ROV
    
    for i_REC = 1:TRIOS_ROV.N_rec(i_ROV)
        
    %if ( zI1f(i_ROV) > 0 ); % ne travailler que sur les données non-nulles
    
        for i_TAG = 1:16;
       
            zaddr_valid = find( zLL(i_TAG,:) > 0 ); % addresses des points valides 
            zdiff       = abs( 	datenum(TRIOS_ROV.date(i_ROV,i_REC))- zDate_TAG(i_TAG,zaddr_valid) ); % différences de dates
            zaddr_min   = find( zdiff < 1e-5); % find difference less than a minute
       
            if ( size(zaddr_min,2) > 0 )
       
                datetime(zDate_TAG(i_TAG,zaddr_min),'ConvertFrom','datenum');
                TRIOS_ROV.date(i_ROV,i_REC);
       
                zaddr_tag_station(i_ROV,i_REC) = i_TAG;
                zaddr_tag_index(i_ROV,i_REC)   = min(zaddr_min);
                zLL_ROV(i_ROV,i_REC) = zLL( i_TAG, zaddr_tag_index(i_ROV,i_REC) ); 
                % on a trouvé la donnée du TAG correspondant Ã  l'échantillonnage ROV
       
            end
        end
    end
end

zaddr_valid = find (TRIOS_ROV.EPAR > 1.0e-5 & ~isnan(TRIOS_ROV.EPAR) & ~isnan(zLL_ROV) );
figure; hold on; box on
plot (zLL_ROV(zaddr_valid),(log10(TRIOS_ROV.EPAR(zaddr_valid))), 'ksq', 'MarkerFaceColor', 'b');
ylabel('log10(EPAR)')
xlabel('LL1+LL2')
% plot (log10(TRIOS_ROV.I1f), zLL1_ROV, 'ksq', 'MarkerFaceColor', 'b');

x = zLL_ROV(zaddr_valid);
y = log10(TRIOS_ROV.EPAR(zaddr_valid));
[p,s] = polyfit (x, y, 1)
plot ([1000:1700], p(2)+p(1)*[1000:1700], 'k:', 'LineWidth', 2)
err = std (y./x) % erreur pente
err1 = std(y) % erreur ordonnée Ã  l'origine
r = corrcoef (x, y) % coefficient de corrélation linéaire
R = r.*r % coefficient de détermination

yfit = polyval(p,x);
resnorm = sum( (yfit-y).^2 );
RMSE    = 10.^(sqrt(resnorm) / size(yfit,1))
% root-mean-square-error des valeurs de EPAR retrouvées à partir de LL1 et LL2


set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);


clear all;
close all;

%---------------------------------------------------------------------------
%
% calibrate_TAG_LARM : test calibration of TAG vs LARM data
%
% start: Martin Vancoppenolle, May 2019
%
%---------------------------------------------------------------------------

user = 'Martin' % 'Bastien'

if ( strcmp(user,'Martin') == 1 )
   indir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/MAT_FILES/'
end

% Ici dessous tu dois changer ton directory...
if ( strcmp(user,'Bastien') == 1 )
   indir = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/MAT_FILES/'
end


% Load tag and trios data
% ATTENTION ICI TU DOIS T'ASSURER QUE LES FICHIERS TAG et LARM sont bien dans "indir",
% sinon �a ira pas
TAG        = load( [ indir, 'PSS117_TAG.mat' ] );
LARM_TRIOS = load( [ indir, 'PSS117_TRIOS_LARM.mat' ] );

% tu dois trouver un truc pour plotter
% LARM_TRIOS.I1f (W/m2) vs TAG.LL1 (dimensionless)
% LARM_TRIOS.I2f (W/m2) vs TAG.LL2 (dimensionless)

% tu ne dois plotter les choses QUE lorsqu'on a une donn�e LARM

% � la fin on veut deux plot 2D avec 
% en x: LARM_TRIOS.I1f et en y: TAG.LL1 
%       LARM_TRIOS.I2f          TAG.LL2

% Les dates LARM_TRIOS.DateNum et TAG.DateNum vont te permettre de faire la correspondance

% un petit exemple ci-dessous te le fait pour zI1f

% ensuite, il faut d�terminer le coefficient "a", et les R2 et tout ce que tu trouves
% d'int�ressant...

zI1f        = reshape( LARM_TRIOS.I1f, 11*4*3, 1 ); % d�plie la matrice en un vecteur 1D
zI2f        = reshape( LARM_TRIOS.I2f, 11*4*3, 1 ); % d�plie la matrice en un vecteur 1D

zDate_LARM  = reshape( LARM_TRIOS.DateNum, 11*4*3, 1 ); % deplie

zLL1        = TAG.LL1;
zLL2        = TAG.LL2;

zDate_TAG   = TAG.DateNum;

zLL1_LARM = nan(11*4*3,1); % le vecteur qui contient les donn�es du TAG dans l'espace des donn�es du TRIOS
zLL2_LARM = nan(11*4*3,1);

for i_LARM = 1:11*4*3 % boucle sur tous les �chantillonnages LARM

    if ( zI1f(i_LARM) > 0 ); % ne travailler que sur les donn�es non-nulles
    
        for i_TAG = 1:16;
       
            zaddr_valid = find( zLL1(i_TAG,:) > 0 ); % addresses des points valides 
            zdiff       = abs( 	zDate_LARM(i_LARM) - zDate_TAG(i_TAG,zaddr_valid) ); % diff�rences de dates
            zaddr_min   = find( zdiff < 1.0e-5 ) % find difference less than a minute
       
            if ( size(zaddr_min,2) > 0 )
       
                datetime(zDate_TAG(i_TAG,zaddr_min),'ConvertFrom','datenum')
                datetime(zDate_LARM(i_LARM),'ConvertFrom','datenum')
       
                zaddr_tag_station(i_LARM) = i_TAG;
                zaddr_tag_index(i_LARM)   = min(zaddr_min);
                zLL1_LARM(i_LARM) = zLL1( i_TAG, zaddr_tag_index(i_LARM) ); %---> voil�, on a trouv� la donn�e du TAG correspondant � l
                                                                            % echantillonnage L-ARM
                zLL2_LARM(i_LARM) = zLL2( i_TAG, zaddr_tag_index(i_LARM) );
       
            end
       
        end
    end
end

figure; 

p1 = plot( log10(zI1f), zLL1_LARM, 'gsq', 'MarkerFaceColor', 'g'); hold on
p2 = plot( log10(zI2f), zLL2_LARM, 'bsq', 'MarkerFaceColor', 'b'); hold on
plot( [-3.5:0.1:0], 60.*[-3.5:0.1:0]+800 , 'k:')

ylabel( 'LL - MK2' )
xlabel( 'LOG10 (Filtered light intensity) - TRIOS (W/m^2)' )

%legend( [p1, p2], { 'LL1', 'LL2' } )
set(gca, 'FontName', 'Helvetica LT Std')
set(gca,'FontSize',12);
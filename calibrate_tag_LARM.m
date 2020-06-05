clear all;
close all;


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
% sinon ça ira pas
TAG        = load( [ indir, 'PSS117_TAG.mat' ] );
LARM_TRIOS = load( [ indir, 'PSS117_TRIOS_LARM.mat' ] );

% tu dois trouver un truc pour plotter
% LARM_TRIOS.I1f (W/m2) vs TAG.LL1
% LARM_TRIOS.I2f (W/m2) vs TAG.LL2

% tu ne dois plotter les choses QUE lorsqu'on a une donnée LARM

% à la fin on veut deux plot 2D avec 
% en x: LARM_TRIOS.I1f et en y: TAG.LL1 
%       LARM_TRIOS.I2f          TAG.LL2

% Les dates LARM_TRIOS.DateNum et TAG.DateNum vont te permettre de faire la correspondance

% un petit exemple ci-dessous te le fait pour zI1f

% ensuite, il faut déterminer le coefficient "a", et les R2 et tout ce que tu trouves
% d'intéressant...

zI1f        = reshape( LARM_TRIOS.I1f, 11*4*3, 1 );
zDate_LARM  = reshape( LARM_TRIOS.DateNum, 11*4*3, 1 );
zLL1        = TAG.LL1;
zDate_TAG   = TAG.DateNum; %can remove add to date once read_TAG is re-run
for i = 1:16;
    for j = 1:5961
       zDate_TAG(i,j)=addtodate(zDate_TAG(i,j),1900,'y');
    end
end

zLL1_LARM = nan(11*4*3);

for i_LARM = 1:11*4*3
    if ( zI1f(i_LARM) > 0 );
    for i_TAG = 1:16;
       i_TAG
       zaddr_valid = find( zLL1(i_TAG,:) > 0 );
       zdiff       = abs( 	zDate_LARM(i_LARM) - zDate_TAG(i_TAG,zaddr_valid) );
       zaddr_min   = find( zdiff < 1.0e-5 ) % find difference less than a minute
       if ( size(zaddr_min,2) > 0 )
       
       datetime(zDate_TAG(i_TAG,zaddr_min),'ConvertFrom','datenum')
       datetime(zDate_LARM(i_LARM),'ConvertFrom','datenum')
       
       zaddr_tag_station(i_LARM) = i_TAG;
       zaddr_tag_index(i_LARM)   = min(zaddr_min);
       zLL1_LARM(i_LARM) = zLL1( i_TAG, zaddr_tag_index(i_LARM) );
       
       end
    end
    end
end

figure; 

plot( log10(zI1f), zLL1_LARM, 'ksq', 'MarkerFaceColor', 'k')

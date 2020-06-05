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

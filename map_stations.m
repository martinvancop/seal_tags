clear all, close all

%----------------------------------------------------------------------------
% Script test pour realiser une carte des stations (en utilisant une boucle)
%----------------------------------------------------------------------------
%
% MartinV, 7 mai 2020, corrige par Bastien le 11 mai
%

figure('Color','w'); 
hold on;

% --- Realisation de la carte ---

coast = load('coast'); % structure qui contient les cotes

% liste des projections ici:
% https://fr.mathworks.com/help/map/summary-and-guide-to-projections.html
% projection polaire stereo
% axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60]) % definit les axes de la carte

% projection stereographique de Gall (localisee pour une region)
axesm('gstereo','MapLonLimit', [-65, 5],'MapLatLimit',[-80 -60]) % definit les axes de la carte
  
framem on % affiche le cadre
axis   off  % affiche les axes
gridm  on  % affiche la grille
mlabel on % affiche les etiquettes des meridiens
plabel on % affiche les etiquettes des paralleles

geoshow(coast.lat,coast.long, 'DisplayType','polygon', 'FaceColor', [135 206 250]/255.) % affiche le continent antarctique

% Ajout des donnees des stations (latitude / longitude)
N_sta = 18
lat_sta = [-69.06889, -69.1024863, -70.1951665, -70.28556, -70.2991597, -70.33583, -70.3516683, -70.38528, -70.50056, -70.50111, -70.63361, -70.4590423, -71.02668, -71.2415445, -67.6691357, -67.5325338, -65.20111, -65.5224318 ] % données latitude (ordre chronologique)
lon_sta = [-0.09944444, -0.0379608, -11.248467, -9.985833, -10.0443703, -8.935278, -8.975854, -8.752222, -8.185,  -9.166667, -9.901667, -8.3512282, -13.4930235, -19.6462462 ,-46.3136623, -46.3817188, -45.935, -46.113281 ] % données longitude (ordre chronologique)

sta_type = [ 1, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2 ] % 1 pour L-arm, 2 pour ROV 

is_tag   = [ 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ] % 1 si de la donnee tag est presente, 0 dans le cas contraire

for i_sta = 1:N_sta

    % stations de type ROV, avec du tag (losange rouge plein)
    if ( sta_type(i_sta) == 2 & is_tag(i_sta) == 1  ) % losange rouge
        H(1) = scatterm( lat_sta(i_sta), lon_sta(i_sta), 40, 'd', 'filled', 'MarkerFaceColor', 'r') 
    end
    
    % stations de type ROV, sans tag (losange rouge vide)
    if ( sta_type(i_sta) == 2 & is_tag(i_sta) == 0 )  
        H(2) = scatterm( lat_sta(i_sta), lon_sta(i_sta), 80, 'd', 'r') 
    end
    
    % stations de type L-arm, avec tag (rond noir plein)
    if ( sta_type(i_sta) == 1 & is_tag(i_sta) == 1 )
        H(3) = scatterm( lat_sta(i_sta), lon_sta(i_sta), 40, 'o', 'filled', 'MarkerFaceColor', 'k')
    end 
    
    % stations de type L-arm, sans tag (rond noir vide)
    if ( sta_type(i_sta) == 1 & is_tag(i_sta) == 0 ) 
       H(4) = scatterm( lat_sta(i_sta), lon_sta(i_sta), 80, 'o', 'k')
    end
    
end

% Axe, titre et legende
xlabel( 'Longitude', 'color', 'k') % changer la couleur du titre en abscisse
ylabel( 'Latitude' , 'color', 'k') % changer la couleur du titre en ordonnee
title( ['PS117 stations'], 'FontSize', 16)
%legend (H, {'\fontsize {11}Ice (L-Arm) Station', '\fontsize {11}ROV Station', '\fontsize {11}Les symboles vides correpondent à des stations ne présentant pas de données issues du tag',''}, 'Location', 'southoutside')
legend ([H(1) H(3)], {'L-ARM', 'ROV'})
%legend ('boxoff') % Pas de cadre autour de la legende


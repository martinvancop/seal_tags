clear all, close all

%
%??????????????????????????????????????????????????
% Script test pour réaliser une carte des stations
%??????????????????????????????????????????????????
% MartinV, 7 mai 2020
%

figure('Color','w'); hold on;

%--- réalisation de la carte

coast = load('coast'); % structure qui contient les cotes

% liste des projections ici:
% https://fr.mathworks.com/help/map/summary-and-guide-to-projections.html
% projection polaire stereo
% axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60]) % definit les axes de la carte

% projection stéréographique de Gall (localisée pour une région)
axesm('gstereo','MapLonLimit', [-90, 20.],'MapLatLimit',[-80 -60]) % definit les axes de la carte
  
framem on % affiche le cadre
axis on   % affiche les axes
gridm on  % affiche la grille
mlabel on % affiche les étiquettes des méridiens
plabel on % affiche les étiquettes des parallèles

geoshow(coast.lat,coast.long, 'DisplayType','polygon', 'FaceColor', [135 206 250]/255.) % affiche le continent antarctique

%--- ajout des données des stations

N_sta = 2

lat_sta = [ -69.06889, -65 ] % ici tu dois remplacer par tes donnees lat / lon
lon_sta = [ -0.09944444, -15 ] % ici tu dois remplacer par tes données lat / lon

sta_type = [ 1, 2 ] % 1 pour ROV, 2 pour L-arm % ici faut remplir aussi

is_tag   = [ 1, 0 ] % 1 si on a de la donnee tag, 0 si pas % ici faut remplir aussi

for i_sta = 1:N_sta

    % stations de type 'ROV', avec du tag
    if ( sta_type(i_sta) == 1 & is_tag(i_sta) == 1) 
       scatterm( lat_sta(i_sta), lon_sta(i_sta), 30, 'ksq', 'MarkerFaceColor', 'k')
    end
    
    % stations de type 'ROV', sans tag
    % [ ajouter le code ici ]
    
    % stations de type 'L-arm', avec tag
    % [ ajouter le code ici ]
    
    % stations de type 'L-arm', sans tag
    if ( sta_type(i_sta) == 2 & is_tag(i_sta) == 0) 
       scatterm( lat_sta(i_sta), lon_sta(i_sta), 30, 'bo' )
    end
    
    
    
end


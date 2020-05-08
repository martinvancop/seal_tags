clear all, close all

% Ce script importe les donn�es d'un fichier L-ARM
% Je l'ai g�n�r� automatiquement � partir de la boite "import data"


%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/WEDDELL/PS117_Compilation_Martin/PS117_Ice_Station_Data/PS117_L_Arm_Spectra_Klaus/PS117_St3_H2_WaterA
%
% To extend the code to different selected data or a different text file, generate a
% function instead of a script.

% Auto-generated by MATLAB on 2020/05/07 14:08:40

%% Initialize variables.
filename = '/Users/ioulianikolskaia/Boulot/ENSEIGNEMENT/MY_LOCEAN_MASTERS/2019-2020/BASTIEN_ALGAE/DATA/WEDDELL/PS117_Compilation_Martin/PS117_Ice_Station_Data/PS117_L_Arm_Spectra_Klaus/PS117_St3_H2_WaterA';
delimiter = '\t';
startRow = 46;
endRow = 300;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error
% occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post processing code is
% included. To generate code which works for unimportable data, select unimportable cells
% in a file and regenerate the script.

%% Allocate imported array to column variable names
wavelength = dataArray{:, 1};
irradiance = dataArray{:, 2};


%% Clear temporary variables
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;

% On fait le plot!
plot(wavelength,irradiance,'LineWidth', 5)
xlabel('\lambda (nm)'), ylabel('Irradiance (mW/m^2/nm)')
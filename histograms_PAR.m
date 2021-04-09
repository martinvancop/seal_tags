


%------------------------
% 1) lecture des donnees
%------------------------

% 1.1 lire les données LL1 et LL2 du TAG

% 1.2 lire les données If1, If2, EPAR, QPAR du TRIOS

%------------------------------------------------------------------------------------------
% 2) reconstruire If1, If2, EPAR, QPAR à partir du tag, en utilisant les regressions qu'on
% a derivees
%------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------
% 3) histogramme de fréquence (toutes les stations ROV)
%------------------------------------------------------------------------------------------

figure
% --- If1 (TRIOS & TAG), 
subplot(2,2,1), hold on, box on;
%histogram() - TAG
%histogram() - TRIOS
xlabel('I_f^1 (W/m^2)')
set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
% --- If2 (TRIOS & TAG), 
subplot(2,2,2), hold on, box on,
xlabel('I_f^2 (W/m^2)')
%histogram() - TAG
%histogram() - TRIOS

set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

% --- EPAR (TRIOS & TAG), 
subplot(2,2,3), hold on, box on,
xlabel('E_{PAR} (W/m^2)')
%histogram() - TAG
%histogram() - TRIOS

set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);

% --- QPAR (TRIOS & TAG)
subplot(2,2,4), hold on, box on,
xlabel('Q_{PAR} (µE/m^2/s)')
%histogram() - TAG
%histogram() - TRIOS

set(gca, 'FontName', 'Helvetica LT Std'); set(gca,'FontSize',12);
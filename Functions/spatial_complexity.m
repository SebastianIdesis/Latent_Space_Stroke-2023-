function [C_FC] = spatial_complexity(static_FC,numBin)
%% Functional Complexity
% Description:
% Complexity as calculated by the integral from Gamora-Lopez et al. 2016
% Calculates the integral over static FC and FCD histograms
%
% input: 
% arg1 = static_FC, rec_Pattern and numBin
% 
% Jakub Vohryzek 25/04/2022 - jakub.vohryzek@upf.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating complexity of the static FC spectrum
[numAreas,~] = size(static_FC);
Isubdiag_FC = find(triu(ones(numAreas),1));

data_Hist_FC = histogram(abs(static_FC(Isubdiag_FC)),numBin,'Normalization','Probability');


Cm_FC = 2*((numBin-1)/numBin); % normalisation factor that represents the extreme case in which the p(rij) is a Dirac-delta function δm

C_FC = 1 - sum(abs(data_Hist_FC.Values - 1/numBin))./Cm_FC;

end
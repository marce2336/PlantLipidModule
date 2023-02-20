function ScoresIdx = GetIdsMatchingScoresList(HarmEssentialMets, EssMetsIdx, MatchingMets)
%%-------------------------------------------------------------------------
% The 'GetIdsMatchingScoresList' creates a boolean vector by comparing the
% identifiers of a defined metabolite from the LipidModule with the common
% metabolites identified in the ModelTemplate. The match between
% identifiers is indicated by an entry equal to 1.
%
% USAGE:
%   ScoresIdx = GetIdsMatchingScoresList(HarmEssentialMets, EssMetsIdx, MatchingMets)
%
% INPUTS:   
%   EssMetsIdx:             Index of LipidModule metabolite.
%           
%   MatchingMets:           List with equivalent metabolites encountered in ModelTemplate.
%
%   HarmEssentialMets:      List of essential metabolites with harmonized
%                           IDs necessary to make the LipidModule functional.
%
% OUTPUTS:
%   ScoresIdx:              Boolean vector with the common metabolites
%                           identified in the ModelTemplate after comparing
%                           the available identifiers for a defined 
%                           LipidModule metabolite.
%
%--------------------------------------------------------------------------                       
%% 
% Generate matching scores list:
ScoresList = cell(4,2);
ScoresList(:,1) = {'PubChem'; 'ChEBIID'; 'InChI'; 'LIPID MAPS'};
ScoresVector = zeros(4,size(MatchingMets,1));

% Retrieve identifiers for metabolite from LipidModule:
Count = 0;

for j = 8:11
    Count = Count + 1;
    ScoresList(Count,2) = cellstr(HarmEssentialMets(EssMetsIdx,j));
end

% Compare identifiers for all matching metabolites:
for i = 1:size(MatchingMets,1)
    
    for j = 1:size(ScoresList,1)
        GetIdentifier = ScoresList{j,2};
        GetScore = sum(contains(MatchingMets(i,:),GetIdentifier));
        ScoresVector(j,i) = GetScore;
    end
end

% Check if any of the identifiers are matched with the LipidModule
% metabolite:
ScoresIdx = zeros(size(MatchingMets,1),1);

for i = 1:size(ScoresVector,2)
    CheckScores = sum(ScoresVector(:,i) == 1);
    ScoresIdx(i) = CheckScores;
end

ScoresIdx = ScoresIdx ~= 0;

end
function [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, flag)
%% ========================================================================
% The 'CheckUniquenessID' function verify that the ID of a metabolite
% 'not connected' in the LipidModule are not duplicated in the metabolites
% list of the ModelTemplate. If the ID of a metabolite is duplicated, the
% 'LM' prefix is added.
%
% USAGE:
%   [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected)
%
% INPUTS:   
%       HarmEssentialMets:      Harmonized list of essential metabolites 
%                               necessary to  make the LipidModule functional.
%           
%       HarmLipidModule:        COBRA model structure that contains the 
%                               parameters of the LipidModule, with
%                               harmonized metabolites.
%
%       HarmLipidModuleMets:    List of unique harmonized metabolites from 
%                               the LipidModule that will be integrated into 
%                               the ModelTemplate.
%
%       NotConnectedIdx:        Index of list of metabolites from 
%                               LipidModule that have not equivalent 
%                               metabolite in ModelTemplate. 
%
%       EsMetsNotConnected:     List of metabolites from LipidModule that 
%                               have not equivalent metabolite in 
%                               ModelTemplate.
%
%       LmID:                   ID of LipidModule metabolite to be
%                               checked.
%
%       EssMetsIdx:             Index of metabolite in HarmEssentialMets
%                               list being checked.
%
%       IsMetLMU:               Vector containing indexes of metabolites
%                               sharing identical IDs between models.
%
%       flag                    ('y') - Update list of metabolites that were not 
%                                       possible to connect to the ModelTemplate
%
% OUTPUTS:
%       HarmEssentialMets:      List of essential metabolites with unique IDs.
%           
%       HarmLipidModule:        COBRA model structure that contains the 
%                               parameters of the LipidModule, with unique IDs.
%
%       HarmLipidModuleMets:    List of unique metabolites from the
%                               LipidModule with unique IDs.
%
%       NotConnectedIdx:        Updated list of metabolites from 
%                               LipidModule that have not equivalent 
%                               metabolite in ModelTemplate. 
%
%       EsMetsNotConnected:     Updated list of metabolites from LipidModule  
%                               that have not equivalent metabolite in 
%                               ModelTemplate.
%         
% =========================================================================
%%

if ~exist('flag', 'var')
    flag = 'n';
end

[~, ExistHarmMetx] = CompareArrays(LmID, HarmLipidModuleMets);
[~, ExistLmMet] = CompareArrays(LmID, HarmLipidModule.mets);

%Add prefix to differentiate the metabolite
if ~isempty(IsMetLMU)
    NewMetID = strcat('LM_',LmID);
    HarmLipidModuleMets(ExistHarmMetx,1) = cellstr(NewMetID); 
    HarmEssentialMets(EssMetsIdx,1) = NewMetID;
    HarmLipidModule.mets(ExistLmMet) = cellstr(NewMetID);
end

% Update list of metabolites not-connected when necessary:
switch flag
    case {'y'}
        NotConnectedIdx = NotConnectedIdx + 1;
        EsMetsNotConnected(NotConnectedIdx,1) = cellstr(LmID);

    case {'n'}
        NotConnectedIdx = '';
        EsMetsNotConnected = '';
end

end
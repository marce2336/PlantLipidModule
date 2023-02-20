function [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected)
%% ========================================================================
% The 'HarmonizeMatchingMets' function harmonizes the metabolites names of
% the ModelTemplate. For this purpose the data in 'met.Names' field of the 
% ModelTemplate are replaced by the data of 'met.Names' field of the
% LipidModule. The metabolites IDs in LipidModule are in turn replaced by 
% the corresponding IDs from ModelTemplate
%
% USAGE:
%   [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected)
%
% INPUTS:   
%       ModelTemplate:          COBRA model structure provided by the user,
%                               with harmonized 'metNames' field.
%
%       HarmEssentialMets:      List of essential metabolites with
%                               harmonized IDs necessary to make the 
%                               LipidModule functional.
%           
%       HarmLipidModule:        COBRA model structure that contains the 
%                               parameters of the LipidModule, with
%                               harmonized IDs.
%
%       HarmLipidModuleMets:    List of unique metabolites from the
%                               LipidModule with harmonized IDs that will  
%                               be integrated into the ModelTemplate.
%
%       ConnectedIdx:           Index of list of metabolites from 
%                               LipidModule harmonized. 
%
%       EsMetsConnected:        List of equivalent metabolites found
%                               between the models.
%
%       LmID:                   ID of LipidModule metabolite to be
%                               harmonized.
%
%       IdxMatchedMets:         Index of equivalent metabolite identified
%                               in ModelTemplate.
%
% OUTPUTS:
%       ModelTemplate:          COBRA model structure provided by the user,
%                               with new entry in harmonized 'metNames' field.
%
%       HarmEssentialMets:      List of essential metabolites with new
%                               harmonized ID.
%           
%       HarmLipidModule:        COBRA model structure that contains the 
%                               parameters of the LipidModule, with new
%                               harmonized IDs.
%
%       HarmLipidModuleMets:    List of unique metabolites from the
%                               LipidModule with new harmonized IDs.
%
%       ConnectedIdx:           Updated index of list of metabolites from 
%                               LipidModule harmonized. 
%
%       EsMetsConnected:        Updated list of equivalent metabolites found
%                               between the models. 
%         
% =========================================================================                       
%% 
% Replace metabolite Name information in TemplateModel:
TmID = ModelTemplate.mets(IdxMatchedMets);

ModelTemplate.metNames(IdxMatchedMets) = HarmLipidModule.metNames(strcmp(HarmLipidModule.mets,LmID));

% Harmonize metID in LipidModule by replacing with corresponding name in ModelTemplate
HarmLipidModule.mets(strcmp(HarmLipidModule.mets,LmID)) = TmID;
HarmEssentialMets(strcmp(HarmEssentialMets(:,1),LmID)) = TmID;

%If Met is present in 'HarmLipidModuleMets' list, it must be removed:
MetExist = find(strcmp(HarmLipidModuleMets(:,1),LmID), 1);

if ~isempty(MetExist) 
    HarmLipidModuleMets(MetExist) = [];
end

% Update Connected Metabolites list:
ConnectedIdx = ConnectedIdx + 1;
EsMetsConnected(ConnectedIdx,1) = cellstr(LmID);
EsMetsConnected(ConnectedIdx,2) = cellstr(TmID);

end
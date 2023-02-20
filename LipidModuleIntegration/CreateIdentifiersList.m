function ListIdentifiers = CreateIdentifiersList(HarmEssentialMets, EssMetsIdx)
%% ========================================================================
% The 'CreateIdentifiersList' funtion retrieves the identifiers available
% for the LipidModule metabolite for which an equivalent metabolite is
% searched in the ModelTemplate.
%
% USAGE:
%   ListIdentifiers = CreateIdentifiersList(HarmEssentialMets, EssMetsIdx)
%
% INPUTS:   
%   ListIdentifiers:        list with identifiers available for the
%                           LipidModule metabolite for which an equivalent
%                           metabolite is searched in the ModelTemplate.
%
% OUTPUTS:
%   HarmEssentialMets:      List of essential metabolites with harmonized 
%                           IDs necessary to make the LipidModule functional.
%
%   EssMetsIdx:             Index of LipidModule metabolite                               
%           
% =========================================================================                       
%% 
ListIdentifiers = cell(4,4);
    
% Create list of identifiers IDs used in Databases:
ListIdentifiers(:,1) = {'PubChem'; 'ChEBI'; 'InChI'; 'LIPID MAPS'};

% Create list of field names in model for each identifier: 
ListIdentifiers(:,2) = {'metPubChemID'; 'metChEBIID'; 'metInChIString'; 'metLIPIDMAPSID'};

% Obtain identifiers information for LipidModule metabolite:
for i = 8:11
    GetIdentifier = HarmEssentialMets(EssMetsIdx, i);
    ListIdentifiers(i-7,3) = cellstr(GetIdentifier);
end

ListIdentifiers = ListIdentifiers(~cellfun(@isempty, ListIdentifiers(:,3)),:);

end
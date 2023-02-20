function [IdentifiersInField, IdentifiersInNotes] = CheckMetIdentifiers(model)
%% ========================================================================

ListIdentifiers = cell(8,2);
Presence = zeros(8,2);

% Add identifiers IDs used in Databases:
ListIdentifiers(:,1) = {'ChEBI'; 'InChI'; 'PubChem';...
    'LIPID MAPS'; 'KEGG'; 'Charge'; 'Formula'; 'metNotes'};

% Add field names in model for each identifier: 
ListIdentifiers(:,2) = {'metChEBIID'; 'metInChIString'; 'metPubChemID';...
    'metLIPIDMAPSID'; 'metKEGGID'; 'metCharges'; 'metFormulas'; 'metNotes'};

%ListIdentifiers(:,4) = {'0'}; % Add default value to array

% Check the existence of fields for metabolites identifiers:
for i = 1:size(ListIdentifiers,1)
    Presence(i,1) = double(isfield(model, ListIdentifiers(i,2)));
end

% If the fields for ChEBI, InChIString, PubChem, Lipid Maps or KEGG
% compounds don't exist, verify if they are included in "metNotes" field
% and act accordingly:

if Presence(end,1) == 1
    for i = 1:size(ListIdentifiers,1)
        IdentifierID = ListIdentifiers{i,1};
        
        switch IdentifierID
            case {'ChEBI'}
                FindChEBI = contains(model.metNotes, 'CHEBI:');
                if sum(FindChEBI) > 0
                    Presence(i,2) = 1;
                end
                
            case {'InChI'}
                FindInChI = contains(model.metNotes, 'InChI=1');
                if sum(FindInChI) > 0
                    Presence(i,2) = 1;
                end
                
            case {'LIPID MAPS', 'Charge', 'Formula', 'PubChem', 'KEGG'}
                FindIdentifier = contains(model.metNotes, IdentifierID);
                if sum(FindIdentifier) > 0
                    Presence(i,2) = 1;
                end
        end
    end  
end

IdentifiersInField = ListIdentifiers(Presence(:,1)> 0,:);
IdentifiersInNotes = ListIdentifiers(Presence(:,2)> 0,:);

end
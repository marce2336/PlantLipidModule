function [ListIdentifiers_new, IdxMatchedMets_new] = RetrieveIdentifiersData(ListIdentifiers_old, IdentifiersInfo, ModelTemplate, Flag1, Flag2)
%% ========================================================================
% The 'RetrieveIdentifiersData' function retieves identifiers data for
% common metabolites encountered between a defined LipidModule metabolite
% and the list of metabolites in the ModelTemplate.
%
% USAGE:
%   [ListIdentifiers, IdxMatchedMets] = RetrieveIdentifiersData(ListIdentifiers, ListOtherIdsInNotes, ModelTemplate, Flag2)
%
% INPUTS:   
%   IdentifiersInfo:        Struct with the following fields:
%                           *.ListOtherIdsInField - List of identifiers
%                             contained in dedicated fields of the 
%                             ModelTemplate.
%
%                           *.ListOtherIdsInNotes - List of identifiers
%                             contained in the 'metNotes' field of the
%                             ModelTemplate.
%
%   ListIdentifiers_old:    List with identifiers available for the
%                           LipidModule metabolite for which an equivalent
%                           metabolite is searched in the ModelTemplate.
%
%   ModelTemplate:          COBRA model structure provided by the user,
%                           with harmonized 'metNames' field.
%
%   Flag1                  (1) ChEBI, InChI, PubChem, and/or LIPID MAPS
%                              identifiers are available in dedicated fields.
%                          (0) Default - Identifiers are not available
%                              in dedicated fields.
%
%   Flag2                  (1) The identifiers information is located
%                              in the 'metNotes' field.
%                          (0) Default - Identifiers are not located in
%                              'metNotes' field.
%
% OUTPUTS:
%   IdxMatchedMets_new:     List of common ('connected') metabolites found
%                           between the models.
%
%
%   ListIdentifiers_new:    List with identifiers available for the
%                           LipidModule metabolite with scores after
%                           comparing againts ModelTemplate identifiers.
%
% =========================================================================                       
%% 
% Extract variables from struct:
ListOtherIdsInNotes = IdentifiersInfo.Notes;
ListOtherIdsInField = IdentifiersInfo.Field;

% First check if identifiers are in dedicated fields and find matching
% mets:
if Flag1 == 1
    IsIDPresent = ismember(ListIdentifiers_old(:,2), ListOtherIdsInField);
    IsIDPresent = ListIdentifiers_old(IsIDPresent, :);

    % If the identifier field in ModelTemplate exist for the identifier
    % available for LipidModule metabolite, search for coincidences:
    if ~isempty(IsIDPresent)
        %Create vector for storing matched metabolites Idxs:
        IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(IsIDPresent,1));
        
        for j = 1:size(IsIDPresent,1)
            ModelField_i = IsIDPresent{j,2};
            IdentifierValue_i = IsIDPresent{j,3};
            MatchIdentifier = contains(ModelTemplate.(ModelField_i), IdentifierValue_i);
            IdxMatchedMets_old(:,j) = MatchIdentifier;

            % If a matching identifier is found, record in scores list:
            if sum(MatchIdentifier) > 0
                Idx = strcmp(ListIdentifiers_old(:,2), ModelField_i);
                ListIdentifiers_old{Idx,4} = '1';
            end
        end

        % If no matching metabolites are found, verify if other identifiers are available in 'met.Notes' field:
        if sum(sum(IdxMatchedMets_old)) == 0 && Flag2 == 1
            %Create vector for storing matched metabolites Idxs:
            IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(ListIdentifiers_old,1));
            
            % Search for coincidences in 'met.Notes' field of ModelTemplate:
            [ListIdentifiers_old, IdxMatchedMets_old] = SearchIDsInNotesField(ListIdentifiers_old, ListOtherIdsInNotes, IdxMatchedMets_old);
        end
    
    % If ModelTemplate doesn't have a dedicated field for available identifier of LipidModule 
    % metabolite, search for coincidences in 'met.Notes' if this field exist:     
    else
        if Flag2 == 1
            %Create vector for storing matched metabolites Idxs:
            IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(ListIdentifiers_old,1));
            
            % Search for coincidences in 'met.Notes' field of ModelTemplate:
            [ListIdentifiers_old, IdxMatchedMets_old] = SearchIDsInNotesField(ListIdentifiers_old, ListOtherIdsInNotes, IdxMatchedMets_old);
        else
            IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(IsIDPresent,1));
        end
    end
    
% If identifiers are not in dedicated fields, but 'met.Notes' contains
% identifiers information, act accordingly:    
elseif Flag2 == 1
    %Create vector for storing matched metabolites Idxs:
    IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(ListIdentifiers_old,1));

    % Search for coincidences in 'met.Notes' field of ModelTemplate:
    [ListIdentifiers_old, IdxMatchedMets_old] = SearchIDsInNotesField(ListIdentifiers_old, ListOtherIdsInNotes, IdxMatchedMets_old);

else
    IdxMatchedMets_old = zeros(size(ModelTemplate.mets,1),size(ListIdentifiers_old,1));
end

ListIdentifiers_new = ListIdentifiers_old;
IdxMatchedMets_new  = IdxMatchedMets_old;
end

% The 'SearchIDsInNotesField' function searches for coincidences between
% identifiers available for LipidModule metabolite and the information
% available in 'met.Notes' field of ModelTemplate:

function [ListIdentifiers_new, IdxMatchedMets_new] = SearchIDsInNotesField(ListIdentifiers_old, ListOtherIdsInNotes, IdxMatchedMets_old)

for k = 1:size(ListIdentifiers_old,1)
    
    % Get identifier available for LipidModule metabolite:
    Identifier_i = ListIdentifiers_old(k,3);
    MatchIdentifier = contains(ListOtherIdsInNotes,Identifier_i);
    IdxMatchedMets_old(:,k) = sum(MatchIdentifier,2);

    if sum(sum(MatchIdentifier)) > 0
        ListIdentifiers_old{k,4} = '1';
    end
end

ListIdentifiers_new = ListIdentifiers_old;
IdxMatchedMets_new  = IdxMatchedMets_old;
end
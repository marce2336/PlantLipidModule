function CreatePairedMetsFile(EsMetsNotConnected, EsMetsConnected, IntegLmMets, ModelTemplate, ListOtherIdsInField, ListOtherIdsInNotes, Flag1, Flag2, Flag3)
%% ========================================================================
% The 'CreatePairedMetsFile' function generates .txt files with the lists of common
% metabolites ('connected') between the LipidModule and the ModelTemplate
% and the metabolites from the LipidModule not encountered in the
% ModelTemplate ('New Metabolites Added').
%
% USAGE:
%   CreatePairedMetsFile(EsMetsNotConnected, EsMetsConnected, IntegLmMets, ModelTemplate, ListOtherIdsInField, ListOtherIdsInNotes, Flag1, Flag2, Flag3)
%
% INPUTS:   
%   ModelTemplate:          COBRA model structure provided by the user,
%                               with harmonized 'metNames' field.
%
%   EsMetsConnected:        List of common ('connected') metabolites found
%                           between the models.
%
%   EsMetsNotConnected:     Metabolites from the LipidModule not encountered
%                           in the ModelTemplate ('New Metabolites Added')
% 
%   IntegLmMets:            Struct that contains the following fields:           
%                           * .HarmLmMetsInteg - list with harmonized metabolite IDs that results
%                             from matching the information of the metabolites
%                             of the LipidModule and ModelTemplate.
%           
%                           * .HarmEssMetsInteg - definite list of the harmonized
%                             metabolites that will be added to the ModelTemplate.
% 
%   ListOtherIdsInField:	List of identifiers contained in dedicated 
%                           fields of the ModelTemplate.
%
%   ListOtherIdsInNotes:    List of identifiers contained in the 'metNotes'
%                           field of the ModelTemplate.
%
%   Flag1                   (1) ChEBI, InChI, PubChem, and/or LIPID MAPS
%                               identifiers are available in dedicated fields.
%                           (0) Default - Identifiers are not available
%                               in dedicated fields.
%
%   Flag2                   (1) The identifiers information is located
%                               in the 'metNotes' field.
%                           (0) Default - Identifiers are not located in
%                               'metNotes' field.
%
%   Flag3                   (1) KEGGIDs are available in 'metKEGGID' field.
%                           (0) KEGGIDs are available in 'metNotes' field.
%
% =========================================================================                       
%% 
% Initialize variables:
if ~exist('Flag3', 'var')
    
    % Double check if KEGGID field is in model:
    if isfield(ModelTemplate, 'metKEGGID')
        Flag3 = 1;
    else
        Flag3 = 0;
    end
end

% Create output files:

% 1). Create output file with metabolites that were not encountered in the ModelTemplate and were newly added!
EsMetsNotConnected = EsMetsNotConnected(~cellfun(@isempty, EsMetsNotConnected(:,1)),:);
NotConIdx = ismember(IntegLmMets.HarmEssMetsInteg(:,1),EsMetsNotConnected);
NewMetsAdded = splitvars(table(cellstr(IntegLmMets.HarmEssMetsInteg(NotConIdx,[1:4,7]))));
NewMetsAdded.Properties.VariableNames = {'Abbreviation','Description','Charged_formula','Charge','KEGG_ID'};
pathMetsAdded = fullfile('OutputFiles','NewMetsAdded.txt');
writetable(NewMetsAdded,pathMetsAdded)
clearvars NotConIdx pathMetsAdded NewMetsAdded

% 2). Create output file with common ('connected') metabolites between the ModelTemplate and the LipidModule:
EsMetsConnected = EsMetsConnected(~cellfun(@isempty, EsMetsConnected(:,1)),:);
PairedMetsIdxLm = ismember(IntegLmMets.HarmEssMetsInteg(:,1),EsMetsConnected(:,2));
PairedMetsLm = cellstr([EsMetsConnected(:,1),IntegLmMets.HarmEssMetsInteg(PairedMetsIdxLm,[2:4,7])]);
PairedMetsIdxTm = ismember(ModelTemplate.mets,EsMetsConnected(:,2));
PairedMetsTm = ModelTemplate.mets(PairedMetsIdxTm);
PairedMetsTm = [PairedMetsTm, ModelTemplate.metNames(PairedMetsIdxTm)];

% If there are common ('connected') metabolites, extract the corresponding
% information and create file:
if ~isempty(PairedMetsTm)
    % Add KEGGID if field is available in model:
    if Flag3 == 1
        GetIds = ModelTemplate.metKEGGID(PairedMetsIdxTm);
        PairedMetsTm = [PairedMetsTm, GetIds];
        LmNames = {'Lm_Abbreviation','Lm_Description','Lm_Charged_formula','Lm_Charge',...
        'Lm_KEGG_ID','Tm_Abbreviation','Tm_Description','Tm_KEGG_ID'};
    else
        LmNames = {'Lm_Abbreviation','Lm_Description','Lm_Charged_formula','Lm_Charge',...
        'Lm_KEGG_ID','Tm_Abbreviation','Tm_Description'};
    end

    % Add identifiers information to ModelTemplate metabolites Paired-up:
    % Get identifiers from dedicated fields:
    FieldNames = '';
    if Flag1 == 1
        for i = 1:size(ListOtherIdsInField,2)
            GetIds = ModelTemplate.(ListOtherIdsInField{i})(PairedMetsIdxTm);
            PairedMetsTm = [PairedMetsTm, GetIds];
            FieldNames = ListOtherIdsInField(i);
        end
    end

    % Get identifiers from 'met.Notes' field:
    if Flag2 == 1 && ~isempty(EsMetsConnected)
        GetIds = ListOtherIdsInNotes(PairedMetsIdxTm, :);
        IsField = sum((strlength(GetIds) <= 6),1);
        IsField = IsField == size(GetIds,1);
        GetIds = GetIds(:,IsField == 0);
        PairedMetsTm = [PairedMetsTm, GetIds];
        GetNotesNames = extractBetween(GetIds(1,:), 1, 5);
        FieldNames = [FieldNames, GetNotesNames];
        FieldNames = strrep(FieldNames, ':', '');
        FieldNames = strrep(FieldNames, '=', '');
    end

    PairedMets = strings(size(EsMetsConnected,1), size(PairedMetsLm,2) + size(PairedMetsTm,2));
    CountPaired = 0;

    for PairedIdx = 1:size(EsMetsConnected,1)
        CountPaired = CountPaired + 1;
        PairedID = EsMetsConnected{PairedIdx,2};
        PairedTm = strcmp(PairedMetsTm,PairedID);
        PairedMets(CountPaired,:) = [(PairedMetsLm(PairedIdx,:)),(PairedMetsTm(PairedTm,:))]; 
    end

    FieldNames = [LmNames, FieldNames];
    PairedMets = splitvars(table(PairedMets));
    PairedMets.Properties.VariableNames = FieldNames;
    pathPairedMets = fullfile('OutputFiles','PairedMets.txt');
    writetable(PairedMets,pathPairedMets)
    clearvars PairedIdx PairedID PairedTm pathPairedMets PairedMetsLm PairedMetsTm PairedMets

% If there are no common ('connected') metabolites, print message:
else
    pathPairedMets = fullfile('OutputFiles','PairedMets.txt');
    fID = fopen(pathPairedMets,'wt');
    fprintf(fID,...
        'No common metabolites were found between the models.\nIt is recommended to add additional information for metabolites, e.g., charge, molecular formula, and KEGGID. \nOther identifiers that are suggested to be added are ChEBI, PubCHEM, InChI string and / or LIPID MAPS');
    fclose(fID);
end
end
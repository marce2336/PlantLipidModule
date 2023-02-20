function [IntegLmMets, HarmLipidModuleInteg, ModelTemplate2] = FindMatchingMets(ModelTemplate, HarmLmVars, IdentifiersInfo, PairedCompartments, UnpairedCompartments)
%% ========================================================================
% The "FindMatchingMets" function search for equivalent metabolites in the
% LipidModule and the ModelTemplate, by comparing the KEGGID, and
% metabolites Charge and Formula, besides other identifiers available. 
% When an identical metabolite is detected, the metabolite ID of the 
% LipidModule is replaced by the corresponding ID of the ModelTemplate.
%
% USAGE:
%   ListMatchedMets = FindMatchingMets(ModelTemplate, HarmEssentialMets, ListOtherIdentifiers, Flag)
%
% INPUTS:   
%        ModelTemplate:         COBRA model structure provided by the user.
%
%        HarmEssentialMets:     List of essential metabolites necessary to
%                               make the LipidModule functional.
%           
%        HarmLipidModule:       COBRA model strcuture that contains the 
%                               parameters of the LipidModule.
%
%        HarmLipidModuleMets:   List of unique metabolites from the
%                               LipidModule that will be integrated into the
%                               ModelTemplate.
%
%        IdentifiersInfo:       Struct with the following fields:
%
%               * .ListOtherIdsInField - List of identifiers contained in 
%                   dedicated fields of the ModelTemplate.
%
%               * .ListOthers - List of identifiers contained in the
%                   'metNotes' field of the ModelTemplate.
%
%               * .Flag1 - (1) ChEBI, InChI, PubChem, and/or LIPID MAPS
%                              identifiers are available in dedicated fields.
%                          (0) Default - Identifiers are not available
%                              in dedicated fields.
%
%               * .Flag2 - (1) The identifiers information is located
%                              in the 'metNotes' field.
%                          (0) Default - Identifiers are not located in
%                              'metNotes' field.
%
%               * .Flag3 - (1) KEGGIDs are available in 'metKEGGID' field.
%                          (0) KEGGIDs are available in 'metNotes' field.
%
%        PairedCompartments:    List of compartments from LipidModule
%                               paired-up with the ModelTemplate compartments.
%
%        UnpairedCompartments:  List of compartments of LipidModule which
%                               don't have a equivalent compartment in
%                               ModelTemplate.
%                              
% OUTPUTS:
%       IntegLmMets:      Struct that contains the following fields:           
%                         * .HarmLmMetsInteg - list with harmonized metabolite IDs that results
%                             from matching the information of the metabolites
%                             of the LipidModule and ModelTemplate.
%           
%                         * .HarmEssMetsInteg - definite list of the harmonized
%                             metabolites that will be added to the ModelTemplate.
%
%        HarmLipidModuleInteg:    Lipid module with harmonized metabolites
%                                 IDs.
%           
%        ModelTemplate2:          COBRA model strcuture of ModelTemplate
%                                 with harmonized 'metNames' field.
% =========================================================================                       
%% 
% Extract variables from structure fields:
HarmLipidModuleMets = HarmLmVars.HarmLmMets; 
HarmEssentialMets   = HarmLmVars.HarmEssMets; 
HarmLipidModule     = HarmLmVars.HarmLM;
ListOtherIdsInNotes = IdentifiersInfo.Notes;
ListOtherIdsInField = IdentifiersInfo.Field;

% Initialize variables:
if ~isfield(IdentifiersInfo, 'Flag1')
    Flag1 = 0;
else
    Flag1 = IdentifiersInfo.Flag1;
end

if ~isfield(IdentifiersInfo, 'Flag2')
    Flag2 = 0;
else
    Flag2 = IdentifiersInfo.Flag2;
end

if ~isfield(IdentifiersInfo, 'Flag3')
    Flag3 = 0;
else
    Flag3 = IdentifiersInfo.Flag3;
end

ListLipidClasses = {'[GL' '[SP' '[GP'};
ConnectedIdx = 0;
NotConnectedIdx = 0;
EsMetsNotConnected = cell(size(HarmEssentialMets,1),1);
EsMetsConnected = cell(size(HarmEssentialMets,1),2);
    
for EssMetsIdx = 1:size(HarmEssentialMets,1)
    
    % Obtain information for LipidModule metabolite:
    LmKEGGidx     = HarmEssentialMets(EssMetsIdx,7); % Met KEGGID
    LmNeutralForm = HarmEssentialMets(EssMetsIdx,6); % Met Neutral Formula
    LmMetName     = HarmEssentialMets(EssMetsIdx,2); % Met Name
    LmID          = HarmEssentialMets(EssMetsIdx,1); % Met ID
    [~, IsMetLMU] = CompareArrays(LmID, ModelTemplate.mets);
    
    % Obtain compartment of LipidModule metabolite:
    LmComp = HarmEssentialMets(EssMetsIdx,5);
    
    % Find equivalente metabolites in ModelTemplate by comparing KEGGIDs:
    if Flag3 == 1
        KEGGmatch = strcmp(ModelTemplate.metKEGGID,LmKEGGidx);
    else
        KEGGmatch = contains(ListOtherIdsInNotes,LmKEGGidx);
        KEGGmatch = sum(KEGGmatch,2);
    end
    
    %% CASE 1: KEGGID is available for LipidModule metabolite, thus it will
    %  be used to find equivalent metabolites in ModelTemplate
    
    if sum(KEGGmatch) > 0 && LmKEGGidx ~= ""
        
        %Build list with matching metabolites from ModelTemplate:
        MatchingMets = string(ModelTemplate.mets(KEGGmatch > 0));
        IdxMatchedMets = find(KEGGmatch);
        
        % Add metabolites Charge and Charged- and Neutral-formulas:
        Count = size(MatchingMets,3);
        MatchingMets = [MatchingMets, strings(size(MatchingMets,1),3)];
        MatchingMets(:,Count+1) = string(ModelTemplate.metCharges(KEGGmatch > 0));
        MatchingMets(:,Count+2) = string(ModelTemplate.metFormulas(KEGGmatch > 0));
        
        % Add data if identifiers are available in dedicated fields:
        if Flag1 == 1
            Count = size(MatchingMets,2);
            MatchingMets = [MatchingMets, strings(size(MatchingMets,1),size(ListOtherIdsInField,1))];
                
            for i = 1:size(ListOtherIdsInField,1)
                GetIds = ModelTemplate.(ListOtherIdsInField{i})(KEGGmatch > 0);
                MatchingMets(:,Count+i) = string(GetIds);
            end
        end
        
        if Flag2 == 1
            Count = size(MatchingMets,2);
            MatchingMets = [MatchingMets, strings(size(MatchingMets,1),size(ListOtherIdsInNotes,2))];

            for i = 1:size(ListOtherIdsInNotes,2)
                GetIds = ListOtherIdsInNotes(KEGGmatch > 0, i);%(IdxList > 0, i);
                MatchingMets(:,Count+i) = string(GetIds);
            end
        end
        
        %Find metabolite location (compartment):
        ListCompartments = GetMetLocation(LmComp, PairedCompartments, UnpairedCompartments);
        
        % Eliminate metabolites from different compartments:
        MatchComp = contains(MatchingMets(:,1), ListCompartments);
        MatchingMets = MatchingMets(MatchComp, :);
        IdxMatchedMets = IdxMatchedMets(MatchComp);

        % Obtain neutral formulas for ModelTemplate metabolites:
        MatchingMets = ObtainMetCharge(MatchingMets, LmNeutralForm);
        
        % Find metabolite(s) with matching neutral formulas:
        if ~isempty(MatchingMets)
            MatchNForms = strcmp(MatchingMets(:,4), LmNeutralForm);
            MatchingMets = MatchingMets(MatchNForms, :);
            IdxMatchedMets = IdxMatchedMets(MatchNForms);
        end
        
        % If the list contains matching metabolite(s), verify identifiers
        % for selected lipid classes:
        IsLipidClass = contains(LmMetName, ListLipidClasses);
        
        if sum(IsLipidClass) > 0
            % Generate matching scores list:
            ScoresIdx = GetIdsMatchingScoresList(HarmEssentialMets, EssMetsIdx, MatchingMets);
            
            % If a matching metabolite is found act accordingly:
            if sum(ScoresIdx) > 0 && ~isempty(MatchingMets)
                % Retrieve metabolite with matching identifiers:
                IdxMatchedMets = IdxMatchedMets(ScoresIdx);
                
                % Harmonize models fields:
                [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate,...
                    ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID,...
                    HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected);
                
            % If no identifier is available the mets cannot be harmonized!.
            % Verify if ModelTemplate contains mets with identical name as 
            % LipidModule met and change if necessary:
            else
                [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                    EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                    HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
            end
        
        % A matching metabolite was found and doesn't belong to selected lipid classes:
        elseif sum(IsLipidClass) <= 0 && ~isempty(MatchingMets)
            
            % Verify if a unique matching metabolite was found:
            
            % When several metabolites in ModelTemplate share same KEGG, 
            % charge and formula, compare other identifiers if available:
            if size(MatchingMets,1) > 1
                ScoresIdx = GetIdsMatchingScoresList(HarmEssentialMets, EssMetsIdx, MatchingMets);
                MatchingMets = MatchingMets(ScoresIdx == 1, :);
                
                % One matching metabolite is retrieved after comparing identifiers information:
                if size(MatchingMets, 1) == 1
                    IdxMatchedMets = IdxMatchedMets(ScoresIdx);
                    [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate,...
                        ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID,...
                        HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected);
                
                % If none or several matching metabolites are still retrieved,
                % the information is not suficient to find the unique equivalente 
                % metabolite in ModelTemplate:
                else
                    [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                        EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                        HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
                end
            
            % Only one metabolite in ModelTemplate share the same KEGG, charge and formula:   
            else
                [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate,...
                    ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID,...
                    HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected);
            end
        else
            [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
        end
    
    %% CASE 2: A matching metabolite was not found after comparing KEGGs and/or the
    %  metabolite in LipidModule has not KEGG assigned. Compare mets IDs and
    %  act accordingly:
    
    % Here it must be checked if other IDs are available besides KEGGID!
    else
        % Get identifiers data for LipidModule metabolite:
        ListIdentifiers = CreateIdentifiersList(HarmEssentialMets, EssMetsIdx);
        
        % If at least one identifier is available for LipidModule
        % metabolite, search for matching identifiers in ModelTemplate:
        if sum(sum(~cellfun(@isempty,ListIdentifiers))) > 0
            % Retrieve information from existent identifiers fields:
            [ListIdentifiers, IdxMatchedMets] = RetrieveIdentifiersData(ListIdentifiers, IdentifiersInfo, ModelTemplate, Flag1, Flag2);
            
            % Check if any of the identifiers were matched and act accordingly:
            IsIDmatched = strcmp(ListIdentifiers(:,4), '1');
            
            % The identifiers are matched!
            if sum(IsIDmatched) > 0
                IdxMatchedMets = find(sum((IdxMatchedMets),2));
                MatchingMets = strings(size(IdxMatchedMets,1), 4);
                MatchingMets(:,1) = ModelTemplate.mets(IdxMatchedMets);
                MatchingMets(:,2) = string(ModelTemplate.metCharges(IdxMatchedMets));
                MatchingMets(:,3) = string(ModelTemplate.metFormulas(IdxMatchedMets));
                
                % Obtain neutral formulas for ModelTemplate metabolites:
                MatchingMets = ObtainMetCharge(MatchingMets, LmNeutralForm);
                
                % Find metabolite(s) with matching neutral formulas:
                MatchNForms = strcmp(MatchingMets(:,4), LmNeutralForm);
                IdxMatchedMets = IdxMatchedMets(MatchNForms);
                
                % Compare the metabolites location:
                %Find metabolite location (compartment):
                ListCompartments = GetMetLocation(LmComp, PairedCompartments, UnpairedCompartments);

                % Eliminate metabolites from different compartments:
                MatchComp = contains(MatchingMets(:,1), ListCompartments);
                MatchingMets = MatchingMets(MatchComp, :);
                
                if ~isempty(IdxMatchedMets)
                    IdxMatchedMets = IdxMatchedMets(MatchComp);
                end
                
                % If a matching metabolite is found act accordingly:
                if ~isempty(IdxMatchedMets)
                    % Harmonize models fields:
                    [HarmLipidModule, HarmEssentialMets, HarmLipidModuleMets, ModelTemplate,...
                        ConnectedIdx, EsMetsConnected] = HarmonizeMatchingMets(HarmLipidModule, LmID,...
                        HarmEssentialMets, HarmLipidModuleMets, ModelTemplate, IdxMatchedMets, ConnectedIdx, EsMetsConnected);

                    % If no identifier is available the mets cannot be harmonized!.
                    % Verify if ModelTemplate contains mets with identical name as 
                    % LipidModule met and change if necessary:
                else
                    [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                        EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                        HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
                end
                
            %==============================================================
            % KEGGID is not available for LipidModule metabolite and the identifiers don't
            % coincide with any metabolite in ModelTemplate!
            else
                [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                        EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                        HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
            end
        %==================================================================
        % The LipidModule metabolite does not have KEGGID or any other identifier available!
        else
            [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                        EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                        HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
        end
    end            
end

HarmLmMetsInteg = HarmLipidModuleMets;
HarmEssentialMetsInteg = HarmEssentialMets;
HarmLipidModuleInteg = HarmLipidModule;
ModelTemplate2 = ModelTemplate;

% Create struct with information of essential metabolites integrated in
% ModelTemplate:
IntegLmMets.HarmEssMetsInteg  = HarmEssentialMetsInteg;
IntegLmMets.HarmLmMetsInteg   = HarmLmMetsInteg;

%% Create output files:

CreatePairedMetsFile(EsMetsNotConnected, EsMetsConnected, IntegLmMets, ModelTemplate, ListOtherIdsInField, ListOtherIdsInNotes, Flag1, Flag2, Flag3)

end
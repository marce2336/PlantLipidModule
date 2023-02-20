function [IntegLmMets,HarmLipidModuleInteg,ModelTemplate2] = FindMatchingIdentifiers(ModelTemplate,...
    HarmLmVars, IdentifiersInfo, PairedCompartments, UnpairedCompartments)
%% ========================================================================
% The "FindMatchingIdentifiers" function search for equivalent metabolites
% in the LipidModule and the ModelTemplate, by comparing the available
% metabolite identifiers information. When an identical metabolite is 
% detected, the metabolite ID of the LipidModule is replaced by the 
% corresponding ID of the ModelTemplate.
%
% USAGE:
%   ListMatchedMets = FindMatchingMets(ModelTemplate, HarmEssentialMets, ListOtherIdentifiers, Flag)
%
% INPUTS: 
%         ModelTemplate:        COBRA model structure provided by the user.
%
%         HarmEssentialMets:    List of essential metabolites necessary to
%                               make the LipidModule functional.
%           
%         HarmLipidModule:      COBRA model strcuture that contains the 
%                               parameters of the LipidModule.
%
%         HarmLipidModuleMets:  List of unique metabolites from the
%                               LipidModule that will be integrated into the
%                               ModelTemplate.
%
%         IdentifiersInfo:      Struct with the following fields:
%
%               * .ListOtherIdsInNotes - List of identifiers contained in the
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
%               * .Flag4 - (1) Formula and Charge parameters are in model fields.
%                          (0) Formula and Charge are not available.
%
%         PairedCompartments:   List of compartments from LipidModule
%                               paired-up with the ModelTemplate
%                               compartments.
%
%         UnpairedCompartments: List of compartments of LipidModule which
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
%       HarmLipidModuleInteg:   Lipid module with harmonized metabolites IDs.
%           
%       ModelTemplate2:         COBRA model strcuture of ModelTemplate
%                               with harmonized 'metNames' field.
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

if ~isfield(IdentifiersInfo, 'Flag4')
    Flag4 = 0;
else
    Flag4 = IdentifiersInfo.Flag4;
end

ConnectedIdx = 0;
NotConnectedIdx = 0;
EsMetsNotConnected = cell(size(HarmEssentialMets,1),1);
EsMetsConnected = cell(size(HarmEssentialMets,1),2);

for EssMetsIdx = 1:size(HarmEssentialMets,1)
    
    % Obtain information for LipidModule metabolite:
    LmID          = HarmEssentialMets(EssMetsIdx,1); % Met ID
    LmNeutralForm = HarmEssentialMets(EssMetsIdx,6); % Met Neutral Formula
    [~, IsMetLMU] = CompareArrays(LmID, ModelTemplate.mets);
    
    % Obtain compartment of LipidModule metabolite:
    LmComp = HarmEssentialMets(EssMetsIdx,5);
     
    % Get identifiers data for LipidModule metabolite:
    ListIdentifiers = CreateIdentifiersList(HarmEssentialMets, EssMetsIdx);
    
    % Search for matching identifiers in ModelTemplate:
    if sum(sum(~cellfun(@isempty,ListIdentifiers))) > 0
        
        % Retrieve information from existent identifiers fields:
        %if Flag1 == 1
        [ListIdentifiers, IdxMatchedMets] = RetrieveIdentifiersData(ListIdentifiers, IdentifiersInfo, ModelTemplate, Flag1, Flag2);
        
        % Check if any of the identifiers were matched and act accordingly:
        IsIDmatched = strcmp(ListIdentifiers(:,4), '1');
        
        % Check that other IDs besides PubChem were matched:
        if sum(IsIDmatched) == 1
            IsPubChem = ListIdentifiers(strcmp(ListIdentifiers(:,1),'PubChem'), 4);
            % It is necessary additional identifiers, besides PubChem to
            % make sure the right equivalent metabolite is selected:
            if strcmp(IsPubChem, '1')
                IsIDmatched = false;
            end
        end
        
        %==================================================================
        % There are metabolites in ModelTemplate with the same identifiers:
        if sum(IsIDmatched) > 0
            IdxMatchedMets = find(sum((IdxMatchedMets),2));
            MatchingMets = ModelTemplate.mets(IdxMatchedMets);
            
            % Compare metabolites location:
            %Find metabolite location (compartment):
            ListCompartments = GetMetLocation(LmComp, PairedCompartments, UnpairedCompartments);
        
            % Eliminate metabolites from different compartments:
            MatchComp = contains(MatchingMets(:,1), ListCompartments);
            MatchingMets = MatchingMets(MatchComp, :);
            IdxMatchedMets = IdxMatchedMets(MatchComp);
            
            %==============================================================
            % Metabolites with the same identifiers and location
            % (compartment) were matched!
            if ~isempty(MatchingMets)
                % If a matching metabolite is found, get corresponding Charge and Formula:
                MatchingMets = [MatchingMets, cell(size(MatchingMets,1),3)];
                
                % The metabolite Charge and Formula is available, thus they 
                % can be used to compare the metabolites:
                if Flag4 == 1
                    MatchingMets(:,2) = cellstr(string(ModelTemplate.metCharges(IdxMatchedMets)));
                    MatchingMets(:,3) = ModelTemplate.metFormulas(IdxMatchedMets);

                    % Obtain neutral formulas for ModelTemplate metabolites:
                    MatchingMets = ObtainMetCharge(MatchingMets, LmNeutralForm);

                    % Find metabolite(s) with matching neutral formulas:
                    MatchNForms = strcmp(MatchingMets(:,4), LmNeutralForm);
                    %MatchingMets = MatchingMets(MatchNForms, :);
                    IdxMatchedMets = IdxMatchedMets(MatchNForms);
                end
                
                %==========================================================
                % If a matching metabolite is found act accordingly:
                if ~isempty(IdxMatchedMets) && size(IdxMatchedMets,1) == 1
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
            % The metabolites location (compartment) is different!
            else
                [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                    EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                    HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
            end
            
        %==================================================================
        % The identifiers available for the metabolite from the LipidModule
        % were not matched to identifiers of TemplateModel!
        else
            [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
                EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
                HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
        end
        
    % There are no identifiers available for the metabolite from the LipidModule!
    else
        [HarmLipidModuleMets, HarmEssentialMets, HarmLipidModule, NotConnectedIdx,...
            EsMetsNotConnected] = CheckUniquenessID(LmID, HarmLipidModuleMets,IsMetLMU,...
            HarmLipidModule, HarmEssentialMets, EssMetsIdx, NotConnectedIdx, EsMetsNotConnected, 'y');
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

CreatePairedMetsFile(EsMetsNotConnected, EsMetsConnected, IntegLmMets, ModelTemplate, ListOtherIdsInField, ListOtherIdsInNotes, Flag1, Flag2)

end
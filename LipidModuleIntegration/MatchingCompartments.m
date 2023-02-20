function MatchingCompartments(modelTemplate, lipidModule, essentialMets)
%% ========================================================================
% The "MatchingCompartments" function identifies compartments in both 
% models and their features: common compartments, missing compartments, 
% and unifies the abbreviations of the common compartments. For pairing up
% the compartments that were not matched in this function, a user interface
% is launched were the user can manually pair up the remaining compartments.
%
% INPUTS:   Template model.
%
%           LipidModule
%
%           The list of Essential metabolites required to make the
%           Lipid module functional after the integration. 
%
% OUTPUT:   The list of paired up compartments that subsequently will be
%           used in the CompHarmonization function.
% =========================================================================                       
%%                       
unpairedComp = cell(size(lipidModule.comps,1),2);%pre-allocate space
harmonizedComp = cell(size(lipidModule.comps,1),4);%pre-allocate space
commonComps = cell(size(lipidModule.comps,1),2);%pre-allocate space
unpairedCount = 0;
harmonizedCount = 0;
commonCount = 0;

pathComps = fullfile('InputFiles','Compartments.txt'); %Load Compartments synonyms list
compartments = cellstr(string(readcell(pathComps)));

compListTemplate = {modelTemplate.comps,modelTemplate.compNames}; %Create array with compartments abbreviations and names from template model
compListLipidModule = {lipidModule.comps,lipidModule.compNames}; %Create array with compartments abbreviations and names from LipidModule

for i = 1:size(compListLipidModule{1},1)
    compAbbLm = compListLipidModule{1}(i,1);
    idxAbbLm = strcmp((compartments(:,1)), compAbbLm);
    compNames = cellstr(compartments(idxAbbLm,2:size(compartments,2)));
    compNames((cellfun(@isempty,compNames)) ==1 ) = [];
    
    % Comparison of compartment synonyms to find existence of equivalent compartment in template model. 
    matched = FindMatchingComps(compNames, compListTemplate{1,2});
    %[~, matched] = CStrAinBP(compNames, compListTemplate{1,2}, 'i');
        if isempty(matched) == 1
            unpairedCount = unpairedCount + 1;
            unpairedComp(unpairedCount,1) = compListLipidModule{1}(i,1);%Add abbreviation of compartment not-matched
            unpairedComp(unpairedCount,2) = compListLipidModule{2}(i,1);%Add name of compartment not-matched
        else
            compAbbTemplate = compListTemplate{1}(matched,1);%Extract compartment abbreviation of Template model
            if strcmp(compAbbLm,compAbbTemplate) == 0
                compAbbLm = strcat((insertBefore(compAbbLm,compAbbLm,"[")),']');%Add brackets to compartment abbreviation LipidModule
                compAbbTemplate = strcat((insertBefore(compAbbTemplate,compAbbTemplate,"[")),']');%Add brackets to compartment abbreviation Template model
                lipidModule.mets(:,1) = strrep(lipidModule.mets,compAbbLm,compAbbTemplate);%Replace compartment abbreviation in the metabolites list of lipidModule
                essentialMets(:,1) = strrep(essentialMets(:,1), compAbbLm,compAbbTemplate); %Replace compartment abbreviation in Essential Mets list
                harmonizedCount = harmonizedCount + 1;
                harmonizedComp(harmonizedCount,1) = compAbbLm;
                harmonizedComp(harmonizedCount,2) = compListLipidModule{2}(i,1);
                harmonizedComp(harmonizedCount,3) = compAbbTemplate;
                harmonizedComp(harmonizedCount,4) = compListTemplate{2}(matched,1);
            else
                commonCount = commonCount + 1;
                commonComps(commonCount,1) = compAbbLm;
                commonComps(commonCount,2) = compListLipidModule{2}(i,1); 
            end 
        end      
end

commonComps = commonComps(~cellfun(@isempty, commonComps(:,1)),:);
pathCommonComps = fullfile('OutputFiles','CommonComps.mat');
save(pathCommonComps, 'commonComps');

unpairedComp = unpairedComp(~cellfun(@isempty, unpairedComp(:,1)),:);%List with compartments not found in the template model
unpairedCompTable = cell2table(unpairedComp);
pathUnpComps = fullfile('OutputFiles','UnpairedCompTable.mat');
save(pathUnpComps, 'unpairedCompTable'); %Input file for FindMatchingComp app

harmonizedComp = harmonizedComp(~cellfun(@isempty, harmonizedComp(:,1)),:);
pathHComps = fullfile('OutputFiles','HarmonizedComp.mat');
save(pathHComps,'harmonizedComp')

harmEssentialMets = essentialMets;
pathHEMets = fullfile('OutputFiles','HarmEssentialMets.mat');
save(pathHEMets,'harmEssentialMets') %Input file for CompHarmonization function

harmLipidModule = lipidModule;
pathHLM = fullfile('OutputFiles','HarmLipidModule.mat');
save(pathHLM,'harmLipidModule') %Input file for CompHarmonization function

templateCompList = [(compListTemplate{1,1}),(compListTemplate{1,2})];
templateCompList=cellfun(@string,templateCompList);
templateCompList = cellstr(join(templateCompList,", "));
pathTmCompList = fullfile('OutputFiles','TemplateCompList.mat');
save(pathTmCompList, 'templateCompList') %Input file for FindMatchingComp app  

pathTm = fullfile('InputFiles','ModelTemplate.mat');
save(pathTm, 'modelTemplate')

%% Run app "FindmatchingComp" to manually match the compartments that were
%not "Paired up" in the previous step. This requires previous installation
%of the file (FindMatchingComp.mlappinstall) in Matlab apps.
matlab.apputil.run('FindMatchingComp'); %Open the "FindMatchingComp" application 

end
function LipidModuleIntegration
%% ========================================================================
% "LipidModuleIntegration" identify compartments in both models and their
% features: common compartments, missing compartments, and performs the 
% harmonization of the compartment abbreviations, followed by the
% integration of the lipid module with the Template model.
% The Lipid Module contains the reconstructed pathways for storage and
% structural lipid synthesis and metabolism in Arabidopsis leaves.
%
% USAGE:
%
% 1. Type "LipidModuleIntegration" in the Command Window.
%
% 2. Select the metabolic model that will be used as a template to
% integrate the lipid module.
% 
% 3. Read carefully the message displayed related to the integration of the
% light reactions in the template model.
%
% 4. A user interface will appear displaying on the left hand side of the
% screen the compartments of the lipid module that were not matched with
% the template model. In case there is a match, select the pair of
% compartments and press Paired up, otherwise select none on the dropdown
% menu on the right hand side and press Pair up. Once this is done, press
% the Finish button.
%
% 5. Once the integration is done, the following message will be displayed:
% "The model was created in "OutputFiles" folder!".
%
% INPUTS:
%    Template model:    Standard COBRA model structure whose metabolite
%                       list must contain the corresponding KEGG
%                       identifiers, the metabolites formula and charge.
%
% OUTPUT:               OutputModel: Standard COBRA model
%                       structure containing the Template model and the
%                       "Lipid module" reactions.
%
%                       OutputModelUnique: Standard COBRA model
%                       structure containing the Template model and the
%                       "Lipid module" reactions, after eliminating the
%                       duplicated reactions.
%
% =========================================================================
%%

%% Find current directory and save path:
currentFolder = pwd;
getDir = what('LipidModuleIntegration');
pathDirSave = fullfile(getDir.path, 'OutputFiles', 'currentFolder.mat');
save(pathDirSave, 'currentFolder')

%% Request user to load a model that will be used as Template:
fileName = [];
model = [];
supportedFileExtensions = {'*.xml;*.sbml;*.sto;*.xls;*.xlsx;*.mat'};
[filepath,~,~] = fileparts(mfilename('fullpath'));
while isempty(fileName) == 1 && isfield(model,'mets') == 0 && isfield(model,'rxns') == 0 && isfield(model,'comps') == 0
    disp('Enter a model: ');
    [fileName, pathName] = uigetfile([supportedFileExtensions, {'Model Files'}], 'Please select the model file');
    
    if isempty(model) == 1 %Load Template model
        cd(pathName)
        model = readCbModel(fileName);
        cd(filepath);
        ModelTemplate = model;
        break
    end
end

pathFile = fullfile('InputFiles','EssentialMetabolites.txt'); %Load Essential metabolites list
delimiter = '\t';
startRow = 2;
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
fileID = fopen(pathFile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
EssentialMets = [dataArray{1:end-1}];

pathLM = fullfile('InputFiles','LipidModule.mat'); %Load LipidModule
Temp1 = load(pathLM);
LipidModule = Temp1.LipidModule;
clearvars pathTm pathFile delimiter startRow formatSpec fileID dataArray...
    ans pathLM Temp1

%Harmonize different usage of brackets.
ModelTemplate.rxns = regexprep(ModelTemplate.rxns,'\(','\[');
ModelTemplate.rxns = regexprep(ModelTemplate.rxns,'\)','\]');
ModelTemplate.rxns = regexprep(ModelTemplate.rxns,'EX_','Ex_');
ModelTemplate.rxns = regexprep(ModelTemplate.rxns,'Sink_','Si_');

%Harmonize subSystems field data type:
if isfield(ModelTemplate, 'subSystems')
    
    % Check structure of field subSystems to identify reactions that
    % belong to several metabolic pathways, whose information is stored
    % in separated cells and act accordingly:
    
    for i = 1:size(ModelTemplate.subSystems, 1)
            SubSystem_i = ModelTemplate.subSystems{i};
            
            % Verify variable class and act accordingly:
            IsClass = class(SubSystem_i);
            
            switch IsClass
                case {'char'}
                    SubSystem_i = cellstr(SubSystem_i);
            end
            
            SubSystem_i = strjoin(SubSystem_i(1,:), ', ');
            ModelTemplate.subSystems{i} = SubSystem_i;
    end
    
    ModelTemplate.subSystems = cellstr(string(ModelTemplate.subSystems));

end

pathTm = fullfile('InputFiles','ModelTemplate.mat');
save(pathTm, 'ModelTemplate')

clearvars Len Len2 Len_i

%% Verify the existence of chloroplast compartment in the model:

[~, ChlPaired] = CompareArrays({'plastid','chloroplast','Plastid_0'}, ModelTemplate.compNames(:,1), 'i');


if isempty(ChlPaired)
    prompt = ['The chloroplast compartment was not found in the model!.' newline 'Do you want to add light reactions to the model? Y/N [Y]: '];
    AddLightRxns = input(prompt,'s');
    if isempty(AddLightRxns)
    AddLightRxns = 'Y';
    end
else
    prompt = ['======================================================================='...
        newline 'The chloroplast compartment exist in the model!.'...
        newline 'NADPH and ATP has to be produced in this compartment!.'...
        newline 'If light reactions are not included, select [Y] to add them. Y/N [Y]: '];
    AddLightRxns = input(prompt,'s');
    if isempty(AddLightRxns)
    AddLightRxns = 'Y';
    end
end
pathAddLight = fullfile('OutputFiles','AddLightRxns.mat');
save(pathAddLight,'AddLightRxns') %Input file for CompHarmonization function
clearvars pathAddLight AddLightRxns prompt

%% Ask if user wants to save extended model in .xml format:
prompt = ['======================================================================='...
    newline 'The extended model will be saved in .xls and .mat formats!.'...
    newline 'Saving in .xml format can take a while depending on the size of the model.'...
    newline 'Do you want to save the extended model in .xml format?. Y/N [N]: '];
flagSave = input(prompt, 's');
if isempty(flagSave)
    flagSave = 'Y';
end
pathFlagSave = fullfile('OutputFiles', 'flagSave.mat');
save(pathFlagSave, 'flagSave')
clearvars pathFlagSave prompt

%%
%Run 'MatchingCompartments' function to match compartmet information
MatchingCompartments(ModelTemplate, LipidModule, EssentialMets);
                
end
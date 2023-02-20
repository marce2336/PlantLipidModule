function TmNeutralForm = ObtainNeutralFormula(TmChargedForm, LmNeutralForm, Charge)
%% ========================================================================
% The "ObtainNeutralFormula" function calculates the neutral formula of a
% metabolite
%
% USAGE:
%   TmNeutralForm = ObtainNeutralFormula(TmChargedForm, LmNeutralForm, Charge)
%
% INPUTS:   
%   TmChargedForm:  Charged formula of metabolite from ModelTemplate. 
%           
%   LmNeutralForm:  Neutral formula of equivalent metabolite from the LipidModule
%
%   Charge:         Corresponding charge of metabolite from ModelTemplate.                    
%
% OUTPUTS:
%   TmNeutralForm:  Neutral formula of the ModelTemplate metabolite 
%                   calculated using the Charged Formula and its respective
%                   Charge.
%           
% =========================================================================                       
%% 

Forms = cellstr([TmChargedForm,LmNeutralForm]);
                        
%Deconstruct the charged Template mets formula
[AtomMatrix,Atoms,~] = atomic(Forms);
NonIdenticalAtoms = find(~(AtomMatrix(:,1) == AtomMatrix(:,2)));

%Adjust atom stoichiometry according to charge
for AtmIdx = 1:length(NonIdenticalAtoms)
    AtomMatrix((NonIdenticalAtoms(AtmIdx,1)),1) = AtomMatrix((NonIdenticalAtoms(AtmIdx,1)),1)-Charge;
end

TmNeutralForm = ''; %Reconstruct neutral formula Template mets

for AtmMatxIdx = 1:size(AtomMatrix,1)
    StoichCoeff = AtomMatrix(AtmMatxIdx,1);
    if StoichCoeff == 1
        StoichStr = "";
    else
        StoichStr = string(StoichCoeff);
    end
    
    %Generate neutral formula for Template mets
    CatenateElement = strcat((Atoms{1,AtmMatxIdx}),StoichStr);
    TmNeutralForm = strcat(TmNeutralForm,CatenateElement);
end

end
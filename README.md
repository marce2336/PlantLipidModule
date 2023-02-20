# Plant Lipid Module

## Description

Collection of scripts that carry out the integration of a comprehensive reconstrution of lipid metabolism of Arabidopsis rosette (referred as the Plant Lipid Module - PLM),
into a COBRA reconstruction of any size (named Template model). The result is a COBRA reconstruction extended with the lipid metabolic network, which is saved in
.xls, .mat and .xml formats.

Last update: 2023-01-04

This repository is administered by Sandra M. Correa (Cordoba@mpimp-golm.mpg.de), Bioinformatics Group, University of Potsdam - Germany


## Installation

### Requirements

- A functional Matlab installation (R2018a or higher)

- The COBRA toolbox for MATLAB <https://opencobra.github.io/cobratoolbox/stable/installation.html>

- Gurobi Optimizer for MATLAB (version 8.0.1 is recommended) <http://www.gurobi.com/registration/download-reg>


### Dependencies (Recommended Software)

- The COBRA toolbox: the instructions to install the COBRA toolbox version for MATLAB can be found in <https://opencobra.github.io/cobratoolbox/stable/installation.html>

- Gurobi solver: the solver and corresponding license can be obtained from <http://www.gurobi.com/registration/download-reg>.
  After downloading the license, add the path to MATLAB

The remaining dependencies listed below are included in the 'Dependencies' folder, and must be installed in advance:

- 'FindmatchingComp' MATLAB app: The files can be found in the subfolder 'FindmatchingComp_App'.
  Initialize MATLAB and then click on the 'APPS' tab located on the MATLAB Toolstrip. Select the option 'Install App'. A graphical user
  interface will be displayed. Select the 'FindMatchingComp.mlappinstall' file.  

- 'StoichTools': The files can be found in the subfolder 'StoichTools'.
  Add the folder and subfolders to the MATLAB path.
  The functions have been adapted to retrieve the individual atoms of chemical formulas for the compounds included in the PLM.
  Jeffrey Kantor (2021). Stoichiometry Tools (https://www.mathworks.com/matlabcentral/fileexchange/29774-stoichiometry-tools),
  MATLAB Central File Exchange. Retrieved October 2, 2021.


### Requirements of the Template Model

- The Template model must have a COBRA structure

- If the Template Model already contains lipid-related reactions, it is recommended to remove them before the integration of the PLM

- Metabolites List requirements: in addition to the minimum fields required in a COBRA structure, the KEGG identifers, charge and charged formulas must be included
  when available. It is also recommended to include ChEBI identifiers

- Check the model for consistency with the COBRA Toolbox function: 'results = verifyModel(model,'simpleCheck', true)'
  In case that a problem is detected in any of the fields of the model, carry out the necessary adjustments


### Usage


1. Install all the dependencies listed above (see section 'Dependencies')

2. Add the path with subfolders for 'LipidModuleIntegration' folder

3. Initialize the CobraToolbox

4. Change the cobra solver to 'gurobi' <changeCobraSolver('gurobi')>

5. Type in the Command Window 'LipidModuleIntegration'

6. Select the model in which the PLM wants to be integrated (an example model is included in the 'ExampleModel' folder)

7. Carefully read the messages displayed in the Command Window reagarding recommendations for the incorporation of light reactions into the Template model, and the
   creation of a copy of the extended model in .xml format

8. Once the PLM is integrated into the Template model, the message 'The model was created in "OutputFiles" folder!' is displayed. The extended model can be found
   in the route "../LipidModuleIntegration/OutputFiles", together with the list of new metabolites added, metabolites harmonized, compartments that were paired up,
   reactions whose direcction was adjusted and compartments for which no match was found in the PLM.
   The model named 'OutputModel' corresponds to the Template model extended with the PLM.
   The model named 'OutputModelUnique' is a more refined version of the above 'OutputModel' where the duplicated reactions were removed.


### Example files:

To integrate the PLM into the AraCore model (doi: 10.1104/pp.114.235358), first make sure that the dependencies required for running the 'LipidModuleIntegration'
software were already installed (see steps 1-2, section 'How to use the software'). Then, follow the instructions below:

1. Initialize the CobraToolbox

2. Add the path with subfolders for 'IntegrateAraCoreModel' and 'ExampleModel' folders

3. Type in the Command Window 'IntegrateAraCoreModel', and read the instructions carefully.
  * a). If the 'IntegrateAraCoreModel' script is executed for the first time, a message will be displayed on the screen regarding the adjustment of the AraCore model
    	 previous to the integration of the PLM: "The adjusted AraCore model doesn´t exist!. Please wait until it is created in the following route...". This is followed by the execution of the 'LipidModuleIntegration' script, where the next message will be displayed: "The 'LipidModuleIntegration' function is going to be launched ...".

  * b). During the execution of the 'LipidModuleIntegration' script read carefully the instructions provided, and pay special attention to the following messages:
    * Message 1: "Enter a model: " The user is requested to enter the AraCore model. The adjusted model can be found in the 'ExampleModel' folder, under the name     'AraCore.mat'

    * Message 2: "The chloroplast compartment exist in the model!. NADPH and ATP has to be produced in this compartment!. If light reactions are not included, select [Y] to add them. Y/N [Y]: ".
    In order to use the AraCore model as example file, the user must enter here: Y

    * Message 3: "The extended model will be saved in .xls and .mat formats!. Saving in .xml format can take a while depending on the size of the model.
		Do you want to save the extended model in .xml format?. Y/N [N]: "
		If the user additionally wishes to save the file in .xml format, please enter here: Y

4. The extended model can be found in the route: "../LipidModuleIntegration/OutputFiles".

5. The extended model can be further refined, by typing 'BalanceExtendedAraCore'
	* The AraCore model extended with the PLM will be subjected to several adjustments, including elimination of unused metabolites, and
	  mass- and charge-balancing of unbalanced reactions.

	* The mass- and charge-balanced extended model can be found in the route: '../IntegrateAraCoreModel/OutputFiles'


## Contributors

-  Correa S. (Cordoba@mpimp-golm.mpg.de), Bioinformatics Group, University of Potsdam - Germany.

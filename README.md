%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MLDSP: Machine Learning with Digital Signal Processing for ultrafast, accurate, and scalable genome classfication at all taxonomic levels

******************************************
*            Setup Instructions          *
******************************************
Note: Setup instructions are to help user setting up the MLDSP project on local machine.

1) Place the downloaded folder MLDSP in MATLAB working directory and add it to the path (right click on MLDSP folder and choose option "Add to Path>Folders and SubFolders").
2) Change MATLAB current working directory to ~/MATLAB/MLDSP/ (Just open the mentioned directoey in MATLAB).

To run the software:
1)Open the workspace directory MLDSP.
2)Click "Run" button under menu option "EDITOR"
To run from the command window, simply type "MainScript"(without"") and press ENTER.


******************************************
*            Common Instructions         *
******************************************

1) MainScript.m - main program that user needs to run for the results.
2) The complete list of provided datasets(used in the paper) can be found under DataBase section of Readme file.
3) The default dataSet is 'Primates'. User can change the dataSet at line number 20 in MainScript.m (e.g. change to dataSet = '****'; for data-set ****).
4) Users can provide their own dataSets (change the dataSet name to TestDataSet (for 4a or 4b) as mentioned in instrunction 3 above),  
a)Using .fasta files
Please make subFolders(representing different clusters) in directory MLDSP/DataBase/TestDataSet and place the respective sequences(.fasta files) in the subfolders. 
b)Using NCBI accesion numbers
For each cluster, make a .txt file with comma separated NCBI accesion numbers; place all .txt files(one per cluster) in directory MLDSP/DataBase/TestDataSet. 


******************************************
*    Customization Instructions/FAQs     *
******************************************
How to use length normalization by max/min/meadian/mean length?
Default is length normalization by median length. All length stats are already calculated and 
user has the option to switch the length normalization parameter.
To switch change the value of variable mLen at line number 31 in MainScript.m 


How to switch among different numerical representations?
Default representation is Purine-pyramidine. User can switch to other numeric representation by making a call to respective function at line number 40 in MainScript.m .
Please look at "Program Modules Description" section of this file to know more about available representations.

How to use different distance measure?
Default is Pearson correlation coefficient (recommended). User can switch to other distance measure. e.g. To use Euclidean distance, change 'cor' to 'euc' at line number 62 in MainScript.m  

How to traverse through 3D plot?
User can rotate the 3D plot to explore from different angles by selecting rotation mode from the menu. User can switch to data cursor mode and select any point by left click to display more information about the selected point.

How to run the software on new data-set?
User can provide .fasta files or .txt files with NCBI accession numbers. Please read instruction 4 under Common Instructions for more details.

******************************************
*       New Sequence Classification      *
******************************************
NOTE:
At present, MLDSP uses supervised machine learning. So, to classify a new sequence, user needs to provide a training dataSet. 
The predicted label will be among the given labels in training dataSet. User can use any of the given dataSets for training or 
can provide own dataSet.
1) Change the dataSet at line number 12 in MainScript.m to the training dataSet name.
2) Run a program (MainScript.m) with common instructions.
3) Use newSeqClassify.m to predict label of a new sequence (3a OR 3b). 
You need to update the *** with the accession number of a new sequence in question:
a) Call newSeqClassify from command window (type following command and press enter):
newSeqClassify('***', mLen, disMat, lg, alabels, clusterNames)  
b) Uncomment call to newSeqClassify.m at line number 132 in MainScript.m before running the MainScript.m.



******************************************
*       Program Modules Description      *
******************************************

MainScript.m : main script that calls other modules.
readFasta.m : function used to read .fasta files from local machine.
downloadFasta.m : function to download .fasta files from NCBI database using given accesion numbers.
lengthCalc.m : function to calculate length stats (maax/min/mean/median length among given sequences).

numMappingAT_CG.m : function for paired numeric representation.
numMappingAtomic.m : function for atomic representation.
numMappingCodons.m : function for codon representation.
numMappingDoublet.m : function for nearnest neighbor based doublet representation.
numMappingEIIP.m : function for EIIP representation
numMappingInt.m : function for integer representation. 
numMappingIntN.m : function for integer(other variant) representation.
numMappingJustA.m : function for justA representation.
numMappingJustC.m : function for justC representation.
numMappingJustG.m : function for justG representation.
numMappingJustT.m : function for justT representation.
numMappingPP.m : function for purine-pyramidine representation (default).
numMappingRandom3.m : function for random representation (among real, purine-pyramidine and justA).
numMappingRandom13.m : function for random representation (among all 13 representations).
numMappingReal.m : function for real representation.

myupdatefcn.m : function to display information about a selected point in 3D plot.
classifictionCode.m : function to perform supervised machine learning classification on given dataset. It produces classification accuracy for all 6 classifiers (4 for more than 2000 sequences) used in the study.
checkDimension.m : function to check dimension of confusion matrix. MATLAB removes zero rows and cols by default from confusion matrix. This function reintroduce zero rows/cols in such cases to keep the size of confusion matrix n x n for given n clusters.
newSeqClassify.m : function to classify a new sequence.

distinguishable_colors.m : function to pick n maximally distinct colors for n clusters. This function is written by Tim Holy and is downloaded from the following link:
https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors

FTM (FATHOM toolbox) : The toolbox provided as part of this distribution is used for distance computation. This open source toolbox is available for free to download.
Jones, D. L. 2015. Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis. College of Marine Science, University of South Florida, St. Petersburg, FL, USA. Available from: http://www.marine.usf.edu/user/djones/ 

************************
*       DataBase       *
************************
Provided DataSets:
3classes, Amphibians, Birds-Fish-Mammals, ClassToSubclass(Actinopterygii), Dengue, DomainToKingdom(Eukaryota-noProtists),
DomainToKingdom(Eukaryota), FamilyToGenus(Cyprinnidae), Fungi, Influenza, Insects, KingdomToPhylum(Animalia),
Mammalia, Mammals, OrderToFamily(Cypriniformes), PhylumToSubphylum(Chordata), Plants, Primates, Protists, 
SubclassToSuperorder(Neopterygii), SubfamilyToGenus(Acheilognathinae), SubphylumToClass(Vertebrata),
SuperorderToOrder(Ostariophysi), Vertebrates

NOTE: 
To reproduce exact results(classification accuracies) from the paper, user should use default MATLAB ClassificationLearner application (Statistics and machine Learning Toolbox) instead of implemented classifictionCode.m
The implemented code will produce similar, but not identical results. Built-in application is slightly faster and easier to use for new user. 
Import variable "ATestlg" to this built-in application and select last column as response, make numberOfFolds 10 and train the classifiers to get classification accuracies (or to reproduce results).
https://www.mathworks.com/help/stats/classification-learner-app.html
To compute distance matrix only(faster processing), code for classification and 3d MDS plot can be commented out.


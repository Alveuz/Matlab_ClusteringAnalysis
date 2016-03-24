clear;clc;

% Load mat files with MicroArray Data (MD) and RegulonDB Operon sets
% converted to a more suitable data type.
% -----------------
% load 'jointMdAndOperonGenes.mat'; %Regulon and MD joint data set
runPath             = what;
tempGenesListPath   = [runPath.path '\MatData\jointGenesList.mat'];
load(tempGenesListPath); % joint genes list for MD and Regulon

%Path for concensus results (txt files)
pathConsensus   = [runPath.path '\xlsFiles\concensusAssignations\'];

[clusIdxesArr, clusTGenes] = clUtils_txt2clusters(pathConsensus);

%% Transform clustering results to its matrix form.
concensusResults = [clusIdxesArr clusTGenes];
cl_1Data_Concensus = clUtils_clusters2matrix( concensusResults, jointGenesList );



aaa = 'bye world';
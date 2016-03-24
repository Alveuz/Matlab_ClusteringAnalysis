clear;clc;

% Load mat files with MicroArray Data (MD) and RegulonDB Operon sets
% converted to a more suitable data type.
% -----------------
% load 'jointMdAndOperonGenes.mat'; %Regulon and MD joint data set
load 'jointGenesList.mat'; % joint genes list for MD and Regulon

%Path for concensus results (txt files)
pathConsensus   = 'C:\Users\Alveuz\Dropbox\Articles\zz_Colaborations\001_Mishael_EColiClustering\Results\ConsensusAssignations\clusters_finales_v3_2\';

[clusIdxesArr, clusTGenes] = clUtils_txt2clusters(pathConsensus);

%% Transform clustering results to its matrix form.
concensusResults = [clusIdxesArr clusTGenes];
cl_1Data_Concensus = clUtils_clusters2matrix( concensusResults, jointGenesList );



aaa = 'bye world';
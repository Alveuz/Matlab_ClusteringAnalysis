clear;clc;

% Load mat files with MicroArray Data (MD) and RegulonDB Operon sets
% converted to a more suitable data type.
% -----------------

runPath             = what;
tempGenesListPath   = [runPath.path '\MatData\Background\jointGenesList.mat'];
operonDataPath      = [runPath.path '\MatData\Background\operonData.mat'];
kmeans456DataPath   = [runPath.path '\MatData\kmeansAssignation\clust456Data.mat'];

load(tempGenesListPath); % joint genes list for MD and Regulon.
load(operonDataPath); % RegulonDB E. coli operonData.
load(kmeans456DataPath); % Assignations of genes for kmeans W/ 456 Centroids.

%Path for concensus results (txt files)
pathConsensus   = [runPath.path '\xlsFiles\concensusAssignations\'];

%% Transform clustering results to its matrix form.
[clusIdxesArr, clusTGenes]  = clUtils_txt2clusters(pathConsensus);
concensusResults            = [clusIdxesArr clusTGenes];
cl_1Data_Concensus          = clUtils_clusters2matrix( concensusResults, jointGenesList );

%% Analyze clustering results for kmeans 456 and Consensus Clustering.
cl_kmeansRes{1,1}   = cl_1Data_Alg456{1,2};
cl_kmeansRes{2,1}   = cl_2Data_Alg456{1,2};
cl_kmeansRes{3,1}   = cl_3Data_Alg456{1,2};
cl_kmeansRes{4,1}   = cl_4Data_Alg456{1,2};
cl_kmeansRes{5,1}   = cl_5Data_Alg456{1,2};
cl_kmeansRes{6,1}   = cl_6Data_Alg456{1,2};
cl_kmeansRes{7,1}   = cl_7Data_Alg456{1,2};
cl_kmeansRes{8,1}   = cl_8Data_Alg456{1,2};
cl_kmeansRes{9,1}   = cl_9Data_Alg456{1,2};
cl_kmeansRes{10,1}  = cl_10Data_Alg456{1,2};

cl_Res{1,1}     = cl_1Data_Concensus{1,2};
cl_Res{2,1}     = cl_kmeansRes;

expSize = length(cl_Res);

for i=1: expSize
    
    cl_iCluster         = cl_Res{i,1};
    tempGenesCl_iSize   = size(cl_iCluster,1);
    
    for j=1:tempGenesCl_iSize
        disp([ 'Experiment :', num2str(i),  ', Comparing Cluster: ', num2str(j)]);
        if(j==1)
            [precCli, recCli] = clusT_PrecRecMeasures( ...
                                        cl_iCluster(j,:), operonVector );
        else
            [tempCliPrec, tempCliRec]  = clusT_PrecRecMeasures( ...
                                        cl_iCluster(j,:), operonVector );
            
            precCli = [precCli; tempCliPrec];
            recCli  = [recCli; tempCliRec];
        end
        
    end
    
    [ FMeasure{1,i}, JaccardIdx{1,i}, Goodness{1,i} ] = clusT_ExtrnlMeasures( ...
                                                    precCli, recCli );
    
    cl_iPrecMatrix{1,i} = precCli;
    cl_iRecMatrix{1,i}  = recCli;
    
    clear precCli;
    clear recCli;
        
end

disp('bye horrible and mean world');
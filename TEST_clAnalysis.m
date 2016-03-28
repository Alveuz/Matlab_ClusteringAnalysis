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

operonMatrix = operonVector;

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

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 2
end

for i = 1: expSize

    tempClustering      = cl_Res{i,1};
    if(~iscell(tempClustering))
        trialsPerClustering = 1;
    else
        trialsPerClustering = size(tempClustering, 1);
    end
    
    switch trialsPerClustering
        case 1
        chosenMethod        = 1;%Calculate Precision/Recall 
                            %between clusterings and RegulonDB Operons.
        [precCli, recCli] = clValidation_ExternalMeasures(...
                                    chosenMethod, ...
                                    tempClustering, ...
                                    operonMatrix );
        
        chosenMethod = 6;%Calculate F-Measure, J-IDX, and G-IDX
                            %between clusterings and RegulonDB Operons.
        [ FMeasure{1,i}, JaccardIdx{1,i}, Goodness{1,i} ] = ...
                                        clValidation_ExternalMeasures( ...
                                                        chosenMethod, ...
                                                        precCli, recCli );
                                                    
        cncnssPrecMatrix{1,i} = precCli;
        cncnssRecMatrix{1,i}  = recCli;
        
        clear precCli;
        clear recCli;
        
        case default
        for j=1:trialsPerClustering
            chosenMethod        = 1;%Calculate Precision/Recall 
                            %between clusterings and RegulonDB Operons.
            cl_iCluster         = tempClustering{j,1};
            if(j==1)
            [precCli, recCli] = clValidation_ExternalMeasures( ...
                                        chosenMethod,...
                                        cl_iCluster, operonMatrix );
            else
            [tempCliPrec, tempCliRec]  = clValidation_ExternalMeasures( ...
                                        chosenMethod,...
                                        cl_iCluster, operonMatrix );

            precCli = [precCli; tempCliPrec];
            recCli  = [recCli; tempCliRec];
            end
            
            chosenMethod = 6;%Calculate F-Measure, J-IDX, and G-IDX
                            %between clusterings and RegulonDB Operons.
            [ FMeasure{1,i}, JaccardIdx{1,i}, Goodness{1,i} ] = ...
                                        clValidation_ExternalMeasures( ...
                                                        chosenMethod, ...
                                                        precCli, recCli );
            
            kmeans456PrecMatrix{1,j} = precCli;
            kmeans456RecMatrix{1,j}  = recCli;
            
            clear precCli;
            clear recCli;
            
        end

    end

end

matlabpool close;           % Close the distributed computing

disp('bye horrible and mean world');
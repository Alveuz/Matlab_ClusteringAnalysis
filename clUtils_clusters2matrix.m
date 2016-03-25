function [ clusteringSummary ] = clUtils_clusters2matrix( clusteringResults, featureList )
% CLUTILS_CLUSTERS2MATRIX
% Using the results of a clustering algorithm, and a list of features
% obtains the corresponding matrix form of size M x N
% where M is the number of clusters of the algoritm,
% and N is the number of features in the list.
%   INPUT:
%   - CLUSTERINGRESULTS is a cell matrix of size M x 2
%   where M is the number of clusters. The first column, 
%   contain a string with the id of the cluster. 
%   The second column, contains a feature assigned to the 
%   aforementioned cluster.
%   - GENESLIST is a cell matrix of size N x 1
%   where N is the number of features which are going to 
%   be considered for the matrix form
%
%   OUTPUT:
%   - CLUSTERING SUMMARY a cell vector of size 1 x 2
%   The first element contains the ids of the clustering's clusters
%   The second element, contains a matrix of size M x N
%   where M is the number of clusters of the algoritm,
%   and N is the number of features in the list.
% Function coded by PhD Guillermo Santamaría-Bonfil.
% alveuz@gmail.com
% https://github.com/Alveuz/Matlab_ClusteringAnalysis

    featuresSize   = size(featureList(:,1),1);

    for i=1:featuresSize
        %Obtain feature of the MD&RgDB joint list
        tempFeat        = cellstr(featureList{i,1});
        %Find this feature in the consensus features list, and get the corresponding
        %idx for the cluster as assigned by the consensus clustering algorithm.
        tempIdx         = find(strcmp(tempFeat, clusteringResults(:,2)));
        if(i==1)
            jointClAndClB    = clusteringResults(tempIdx,:);
        else
            jointClAndClB    = [jointClAndClB; clusteringResults(tempIdx,:)];
        end
    end

    clustersIds = unique(jointClAndClB(:,1));
    cl_iVector  = int8(zeros(1, featuresSize));

    for i=1:length(clustersIds)
        tempIdClus          = clustersIds(i,1);
        temp_cl_i_Idxes     = find(strcmp(tempIdClus, jointClAndClB(:,1)));

        if(length(temp_cl_i_Idxes)>1)
            cli_GenesSize   = size(temp_cl_i_Idxes, 1);
            subGenesSet     = jointClAndClB(temp_cl_i_Idxes, 2);

            for j=1:cli_GenesSize
                if(j==1)
                    tempIdxes = find(strcmp(subGenesSet(j), featureList(:,1))~=0);
                else
                    tempIdxes = [tempIdxes; ...
                        find(strcmp(subGenesSet(j), featureList(:,1))~=0)];
                end
            end
        else
            tempIdxes = find(strcmp(tempIdClus, featureList(:,1))~=0);
        end

        if(i==1)
            cl_iVector(1,tempIdxes) = 1;
        else
            tempVector              =  int8(zeros(1, featuresSize));
            tempVector(1,tempIdxes) = 1;
            cl_iVector              = [cl_iVector; tempVector];
        end

    end

    clusteringSummary{1,1} = clustersIds;
    clusteringSummary{1,2} = cl_iVector;


end


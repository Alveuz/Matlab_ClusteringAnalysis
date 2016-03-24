function [ clusteringSummary ] = clUtils_clusters2matrix( clusteringResults, genesList )
% CLUTILS_CLUSTERS2MATRIX
% Using the results of a clustering algorithm, and a list of genes
% obtains the corresponding matrix form of size M x N
% where M is the number of clusters of the algoritm,
% and N is the number of genes in the list.
%   INPUT:
%   - CLUSTERINGRESULTS is a cell matrix of size M x 2
%   where M is the number of clusters. The first column, 
%   contain a string with the id of the cluster. 
%   The second column, contains a gene assigned to the 
%   aforementioned cluster.
%   - GENESLIST is a cell matrix of size N x 1
%   where N is the number of genes which are going to 
%   be considered for the matrix form
%
%   OUTPUT:
%   - CLUSTERING SUMMARY a cell vector of size 1 x 2
%   The first element contains the ids of the clustering's clusters
%   The second element, contains a matrix of size M x N
%   where M is the number of clusters of the algoritm,
%   and N is the number of genes in the list.

    genesFeaturesSize   = size(genesList(:,1),1);

    for i=1:genesFeaturesSize
        %Obtain gene of the MD&RgDB joint list
        tempGene        = cellstr(genesList{i,1});
        %Find this gene in the consensus genes list, and get the corresponding
        %idx for the cluster as assigned by the consensus clustering algorithm.
        tempIdx         = find(strcmp(tempGene, clusteringResults(:,2)));
        if(i==1)
            jointClAndDBGenes    = clusteringResults(tempIdx,:);
        else
            jointClAndDBGenes    = [jointClAndDBGenes; clusteringResults(tempIdx,:)];
        end
    end

    clustersIds = unique(jointClAndDBGenes(:,1));
    cl_iVector  = int8(zeros(1, genesFeaturesSize));

    for i=1:length(clustersIds)
        tempIdClus          = clustersIds(i,1);
        temp_cl_i_Idxes     = find(strcmp(tempIdClus, jointClAndDBGenes(:,1)));

        if(length(temp_cl_i_Idxes)>1)
            cli_GenesSize   = size(temp_cl_i_Idxes, 1);
            subGenesSet     = jointClAndDBGenes(temp_cl_i_Idxes, 2);

            for j=1:cli_GenesSize
                if(j==1)
                    tempIdxes = find(strcmp(subGenesSet(j), genesList(:,1))~=0);
                else
                    tempIdxes = [tempIdxes; ...
                        find(strcmp(subGenesSet(j), genesList(:,1))~=0)];
                end
            end
        else
            tempIdxes = find(strcmp(tempIdClus, genesList(:,1))~=0);
        end

        if(i==1)
            cl_iVector(1,tempIdxes) = 1;
        else
            tempVector              =  int8(zeros(1, genesFeaturesSize));
            tempVector(1,tempIdxes) = 1;
            cl_iVector              = [cl_iVector; tempVector];
        end

    end

    clusteringSummary{1,1} = clustersIds;
    clusteringSummary{1,2} = cl_iVector;


end


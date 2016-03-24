function [clusIdxesArr, clusTGenes] = clUtils_txt2clusters( pathStr )
% CLUTILS_TXT2CLUSTERS 
% Load clusters txt files resulting from a clustering algorithm.
% For instance, consensus clustering results, each txt file per
% cluster.
%   INPUT:
%   - PATHSTR is a string containing the folder 
%   path where txt files are located.
%
%   OUTPUT:
%   - CLUSIDXESARR a cell vector of size N x 1
%   repeated indexes correspond to the same cluster
%   - CLUSTGENES is a cell vector of size N x 1
%   together with CLUSIDXESARR match to the clusters
%   assignation

    wildCard        = '*.txt';
    filesArr        = dir([pathStr wildCard]);
    filesArrSize    = size(filesArr,1);

    for i=1:filesArrSize
        tempClusId      = cellstr(filesArr(i,1).name);
        fileID          = fopen([pathStr tempClusId{1,1}]);
        tempContent     = textscan(fileID, '%s');
        tempContent     = tempContent{1,1};
        fclose(fileID);

        if(i==1)
            auxClusId       = strsplit(tempClusId{1,1}, '.');
            tempClusId      = cellstr(auxClusId{1,1});

            numOfGenes(i,1) = length(tempContent);

            if(length(tempContent)>1)
                repMatNum       = length(tempContent);
                clusIdxesArr    = cellstr(repmat(tempClusId,repMatNum,1));
    %             clusIdxesArr    = tempClusId;            
    %             tempContent     = tempContent';
                clusTGenes      = tempContent;
            else
                clusIdxesArr    = tempClusId;
                clusTGenes      = tempContent;
            end

        else
            auxClusId       = strsplit(tempClusId{:,1}, '.');
            tempClusId      = cellstr(auxClusId{1,1});                        

            numOfGenes(i,1) = length(tempContent);

            if(length(tempContent)>1)
                repMatNum           = length(tempContent);
                tempClusIdxesArr    = cellstr(repmat(tempClusId,repMatNum,1));
                clusIdxesArr        = [clusIdxesArr; tempClusIdxesArr];
    %             tempContent   = tempContent';
    %             tempContent1{1,1}   = tempContent;
                clusTGenes          = [clusTGenes; tempContent];
    %             clusTGenes          = [clusTGenes; tempContent];
            else
                clusIdxesArr    = [clusIdxesArr; tempClusId];
                clusTGenes      = [clusTGenes; tempContent];
            end
        end
    end



end


function [ result, varargout ] = clValidation_ExternalMeasures( ...
                                            chosenMethod, ...
                                            clusteringAssgntns,...
                                            varargin)
% CLVALIDATION_EXTERNALMEASURES
% This function calculates several external validation 
% measures for clustering. 
% CHOSENMETHOD is an integer between 1-5. It allows to 
% chose one of the five metrics:
%   1) Precision and Recall between two clusterings
%   2) Precision and Recall between two clusters.
%   3) F-Measure        between two clusterings.
%   4) Jaccard Index    between two clusterings.
%   5) Goodness Index   between two clusterings.
% 
% Each measure has different types of inputs and outputs.
% In the following, these are described.
%   1) PREC/REC BETWEEN TWO CLUSTERINGS.
%   INPUT:
%   - CLUSTERINGASSGNTNS a matrix of size       P x N
%   where P are the number of clusters,
%   and N are the number of features (genes).
%   - VARARGIN -> CLUSTERING_B a matrix of size Q x N
%   where Q are the number of base clusters,
%   and N are the number of features (genes).
%
%   OUTPUT:
%   - RESULT -> PRECISION a matrix of size P X N containing 
%   the precision between A and B clusterings.
%   Where P are the number of clusters of clustering A,
%   and N are the number of features (genes).
%   - VARARGOUT -> RECALL a matrix of size P X N containing 
%   the recall between A and B clusterings.
%   Where P are the number of clusters of clustering A,
%   and N are the number of features (genes).
%
%   2) PREC/REC BETWEEN TWO CLUSTERS.
%   INPUT:
%   - CLUSTERINGASSGNTNS a vector of size       1 x N
%   N are the number of genes.
%   - VARARGIN -> CLUSTERING_B a vector of size 1 x N
%   a base cluster of N number of features.
%
%   OUTPUT:
%   - RESULT -> PRECISION a vector of size 1 X N containing 
%   the precision between A and B clusters.
%   N are the number of features (genes).
%   - VARARGOUT -> RECALL a vector of size 1 X N containing 
%   the recall between A and B clusters.
%   N are the number of features (genes).
%
% Function coded by PhD Guillermo Santamaría-Bonfil.
% alveuz@gmail.com
% https://github.com/Alveuz/Matlab_ClusteringAnalysis

    switch chosenMethod
%       1) PREC/REC BETWEEN TWO CLUSTERINGS.
        case 1
            cl_iClustering      = clusteringAssgntns;
            cl_jClustering      = varargin{1,1};
            tempGenesCl_iSize   = size(cl_iClustering,1);
            for j=1:tempGenesCl_iSize
                if(j==1)
                    [precCli, recCli] = ...
                                clValidation_PrecRecMeasures( ...
                                cl_iClustering(j,:), cl_jClustering );
                else
                    [tempCliPrec, tempCliRec]  = ...
                                    clValidation_PrecRecMeasures( ...
                                    cl_iClustering(j,:), cl_jClustering );

                    precCli = [precCli; tempCliPrec];
                    recCli  = [recCli; tempCliRec];
                end
            end

            result          = precCli;
            varargout{1,1}  = recCli;

            clear precCli;
            clear recCli;
        
%       2) PREC/REC BETWEEN TWO CLUSTERS.
        case 2
            cl_iClustering         = clusteringAssgntns;
            cl_jClustering        = varargin{1,1};
            
            [result, varargout{1,1}]  = ...
                                    clValidation_PrecRecMeasures( ...
                                    cl_iClustering, cl_jClustering );

            
        case 3
            precMatrix  = clusteringAssgntns;
            recMatrix   = varargin{1,1};
            
            coeffBeta = precMatrix.*recMatrix;
            coeffZeta = precMatrix+recMatrix;
            FMeasure  = (2.*coeffBeta)./coeffZeta;
            
            FMeasure(isnan(FMeasure))       = 0;

        case 4
            precMatrix  = clusteringAssgntns;
            recMatrix   = varargin{1,1};
            
            coeffBeta   = precMatrix.*recMatrix;
            coeffZeta   = precMatrix+recMatrix;
            JaccardIdx  = (coeffBeta)./sqrt(coeffZeta - coeffBeta);
            
            JaccardIdx(isnan(JaccardIdx))   = 0;

        case 5
            precMatrix  = clusteringAssgntns;
            recMatrix   = varargin{1,1};
            
            coeffBeta = precMatrix.*recMatrix;
            coeffZeta = precMatrix+recMatrix;
            Goodness  = 0.5.*(coeffZeta);
            
            Goodness(isnan(Goodness))       = 0;
            
    end

end


% -------------------------------------------------------------------------
% [Ben] 3/28/18
% Takes in the x,y,z coordinates of cortical BBs and returns estimated
% anterior and posterior poles of the cell.
% -------------------------------------------------------------------------

function [pole1, pole2] = findPoles(x,y,z)
x=x*0.125;
y=y*0.125;
z=z*0.3;


numBBs=length(x);

BBsDistanceMatrix=zeros(numBBs,numBBs);
maxDistanceForEachBB=zeros(numBBs,1);
for i=1:numBBs
    for j=1:numBBs
        BBsDistanceMatrix(i,j)=sqrt((x(i,1)-x(j,1))^2+...
            (y(i,1)-y(j,1))^2+(z(i,1)-z(j,1))^2);
    end
    maxDistanceForEachBB(i,1)=max(BBsDistanceMatrix(i,:));
end

% use numBBFindPoles of max distance to determine poles
numBBFindPoles=10;
maxDistanceForEachBB_sorted=sort(maxDistanceForEachBB,'descend');

BBPairWithMaxDistance_i=zeros(numBBFindPoles,1);
BBPairWithMaxDistance_j=zeros(numBBFindPoles,1);

% pg. 18 of Jingyi's paper. Get pairwise distance matrix for all detected
% cortical BBs. Find the 10 cortical BB pairs that produce the greatest
% inter-BB distance. 
for k=1:numBBFindPoles
    for i=1:numBBs
        for j=1:numBBs
            if BBsDistanceMatrix(i,j)==maxDistanceForEachBB_sorted(k)
                BBPairWithMaxDistance_i(k)=i;
                BBPairWithMaxDistance_j(k)=j;
            end
        end
    end
end

BBPairWithMaxDistance=[BBPairWithMaxDistance_i;BBPairWithMaxDistance_j];
% extracts unique values in matrix
BBPairWithMaxDistance=unique(BBPairWithMaxDistance);

idxBBPairWithMaxDistance=zeros(length(BBPairWithMaxDistance),3);
for i=1:length(BBPairWithMaxDistance)
    idxBBPairWithMaxDistance(i,1)=x(BBPairWithMaxDistance(i));
    idxBBPairWithMaxDistance(i,2)=y(BBPairWithMaxDistance(i));
    idxBBPairWithMaxDistance(i,3)=z(BBPairWithMaxDistance(i));
end

% pg. 18 of Jingyi's paper. Apply k-means to basal bodies forming the 10
% greatest inter-BB distances. Anterior and posterior pole defined as the
% centroids of two clusters. 
numOfGroup=2;
[~, groupCenter]=kmeans(idxBBPairWithMaxDistance,numOfGroup);
pole1 = groupCenter(1, :) ./ [0.125 0.125 0.3];
pole2 = groupCenter(2, :) ./ [0.125 0.125 0.3];

end
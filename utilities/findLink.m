% -------------------------------------------------------------------------
% [Ben] 12/07/17
% Returns BBs that are part of a valid anterior-posterior pair, together
% with their corresponding partner (i.e. the other half making up the pair)
% -------------------------------------------------------------------------


function [startPt, endPt] = findLink(x, y, z, antPole, postPole, dist2Ant)
n = length(x);
pts_start = 1:n;   % start: posterior
pts_end = 1:n;    % end: anterior

% add distance from point i to plane defined by (j,antPole,postPole) to
% pairwise distance between point i and j
% min dist between partner and plane(BB maxima, antPole, postPole)

% So dist_forward seems to be a matrix where entry (i,j) stores the sum of
% the distance between BBs i and j (indexed based on the order of
% getBBIdx's output) and the distance between j and the plane defined by
% the 2 poles (findPoles.m) and i.

% This metric calculated in both the forward and reverse directions and BBs
% are assigned a partner only if the forward and reverse directions
% identify the same reciprocal BB pair.
dist_forward = zeros(n, n);
dist_Eu = zeros(n, n);
for i = 1:n
    for j = 1:n
        pt = [x(i), y(i), z(i)];
        pt1 = [x(j), y(j), z(j)];
        dist_Eu(i, j) = distance_pts(pt, pt1);
        dist_forward(i, j) = dist_Eu(i, j) + ...
            distance_pt2plane(pt1, pt, antPole, postPole);
    end
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% % IN PROGRESS (Ben)
% n = length(x);
% dist_forward = zeros(n, n);
% dist_Eu = zeros(n, n);
% for i = 1:n
%     for j = 1:n
%         if i == j
%             % by right this should be 0, but for later convenience when we
%             % have to take the minimum along a row we want to exclude
%             % considering this value of 0
%             dist_forward(i, j) = Inf;
%         else
%             pt1 = [x(i) y(i) z(i)];
%             pt2 = [x(j) y(j) z(j)];
%             dist_Eu(i, j) = distance_pts(pt1, pt2);
%             dist_forward(i, j) = dist_Eu(i, j) + ...
%                 distance_pt2plane(pt2, pt1, antPole, postPole);
%         end
%     end
% end
% 
% % keeps track of BBs that have already been assigned as an anterior partner
% assigned = repelem(false, n); 
% skip = false;
% startPt = [];
% endPt = [];
% for i = 1:n
%     % potential anterior partners 
%     [distances, ant_partners] = sort(dist_forward(i, :), 'ascend');
%     j = 1; % ant_partners(1) is always going to be i
%     ant_partner = ant_partners(j);
%     while assigned(ant_partner) && (dist2Ant(ant_partner) >= dist2Ant(i))
%         if j < n
%             j = j + 1;
%             ant_partner = ant_partners(j);
%         else
%             skip = true; % no potential partners
%             break
%         end
%     end
%     startPt = [startPt i];
%     if (ant_partner == i) || skip
%         endPt = [endPt 0];
%         skip = false;
%     else
%         min_dist = distances(j);
%         % distance between BB and its anterior partner less than 5um
%         if min_dist < 5
%             [~, pos_partner] = min(dist_forward(ant_partner, :));
%             % metric in two directions identified the same BB pair
%             if pos_partner == i
%                 assigned(ant_partner) = true;
%                 endPt = [endPt ant_partner];
%             else
%                 endPt = [endPt 0];
%             end
%         else
%             endPt = [endPt 0]; 
%         end
%     end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% dist_backward should be a matrix where entry (i,j) stores the sum of the
% distance between BBs j and i, and the distance between i and the plane
% defined by the 2 poles and j. 

% dist_backward=zeros(length(pts_end),length(pts_start));
% for i=1:length(pts_end)
%     for j=1:length(pts_start)
%         pt=[x(pts_end(i)),y(pts_end(i)),z(pts_end(i))];
%         pt1=[x(pts_start(j)),y(pts_start(j)),z(pts_start(j))];
%         pt2=antPole;
%         pt3=postPole;
%         dist_backward(i,j)=dist_Eu(i,j) + distance_pt2plane(pt1,pt,pt2,pt3);
%     end
% end 
dist_backward = dist_forward;

% You find the anterior partner of i by finding the j such that
% dist_forward(i,j) is smallest and dist2Ant(j) < dist2Ant(i).
% You find the posterior parter of i by finding the j such that
% dist_backward(i,j) is smallest and dist2Ant(j) > dist2Ant(i).

startPt = [];
endPt = [];
prevLen = n;

startPt_forward = zeros(n ,1);
endPt_forward = zeros(n, 1);
% finds anterior partners using forward direction results
for i = 1:n
    distances = dist_forward(i,:); % length: len(pts_end)
    [~, idx] = sort(distances, 'ascend'); % gives ordering to rearrange in ascending order
    for j = 1:n
        if (idx(j) ~= i) && ...
                (dist2Ant(idx(j)) < dist2Ant(i)) % partner is anterior neighbor
            startPt_forward(i) = i; % why not just startPt_forward(i) = i?
            endPt_forward(i) = idx(j); 
            break
        end
    end
end

startPt_backward=zeros(n, 1);
endPt_backward=zeros(n, 1);
% finds posterior partners using backward direction results
for i = 1:n
    distances = dist_backward(i, :); % length: len(pts_start)
    [~, idx] = sort(distances, 'ascend');
    for j = 1:n
        if (idx(j) ~= i) && ...
                (dist2Ant(idx(j)) > dist2Ant(i)) % partner is posterior neighbor
            startPt_backward(i) = i;
            endPt_backward(i) = idx(j);
            break
        end
    end
end

% So now endPt_forward contains all the anterior neighbors if they exist,
% and endPt_backward contains all the posterior neighbors if they exist.
% Note that startPt_forward and startPt_backward will have a 0 at entry i
% if BB_i does not have an anterior partner and does not have a posterior
% partner respectively.

% At the end of the connection process, the only unconnected BBs that
% remain are those that will never yield reciprocal connections. These BBs
% are connected to the closest unconnected BB only if the connection is
% less than 5um. 
if ~isempty(startPt_forward)
    for i = 1:length(startPt_forward)
        % BB_i has anterior partner
        if startPt_forward(i) ~= 0 && endPt_forward(i) ~= 0
            % find anterior partner of BB_i
            idx = find(startPt_backward == endPt_forward(i));
            % check that posterior partner of bb_i's anterior partner is in
            % fact BB_i
            if endPt_backward(idx) ~= 0 && endPt_backward(idx) == startPt_forward(i)
                % p1 is BB_i and p2 is p1's anterior partner
                p1 = [x(startPt_forward(i)), y(startPt_forward(i)), z(startPt_forward(i))];
                p2 = [x(endPt_forward(i)), y(endPt_forward(i)), z(endPt_forward(i))];
                % distance between BB and its partner < 5um
                if distance_pts(p1, p2) <= 5
                    % startPt and endPt together record pairs of BBs and
                    % their anterior partners.
                    startPt = [startPt, startPt_forward(i)];
                    endPt = [endPt, endPt_forward(i)];
                end
            end
        end
    end
end
%---
% % Of the BBs that are not part of an anterior-posterior pair, identify
% % the BBs that did not have an valid anterior neighbor as
% % well as those that did not have a valid posterior neighbor.
% % Note that setdiff sorts the output.
% pts_start = setdiff(pts_start, startPt);
% pts_end = setdiff(pts_end, endPt);
% 
% % create dist_forward and dist_backward again, this time based only on
% % BBs that have not already been assigned to a valid anterior-posterior
% % pair
% dist_forward = zeros(length(pts_start), length(pts_end));
% dist_Eu = zeros(length(pts_start), length(pts_end));
% for i = 1:length(pts_start)
%     for j = 1:length(pts_end)
%         pt = [x(pts_start(i)),y(pts_start(i)),z(pts_start(i))];
%         pt1 = [x(pts_end(j)),y(pts_end(j)),z(pts_end(j))];
%         dist_Eu(i,j) = distance_pts(pt, pt1);
%         dist_forward(i,j)=dist_Eu(i,j) + ... 
%             distance_pt2plane(pt1, pt, antPole, postPole);
%     end
% end
% 
% dist_backward = zeros(length(pts_end), length(pts_start));
% for i = 1:length(pts_end)
%     for j = 1:length(pts_start)
%         pt = [x(pts_end(i)), y(pts_end(i)), z(pts_end(i))];
%         pt1 = [x(pts_start(j)), y(pts_start(j)), z(pts_start(j))];
%         dist_backward(i,j) = dist_Eu(i,j) + ...
%             distance_pt2plane(pt1, pt, antPole, postPole);
%     end
% end
% currLen = length(pts_end);
%     
% % Repeats above procedure until no further progress can be made - i.e.
% % every BB has been successfully assigned to a pair, or the remaining
% % inter-BB distances are too large, or no new pair was added (meaning that
% % the metric did not identify the same pair when applied in both the
% % forward and backward direction. 
% while ~isempty(pts_start) && min(min(dist_Eu)) <= 5 && currLen~=prevLen
%     prevLen=currLen;
%     startPt_forward=zeros(length(pts_start),1);
%     endPt_forward=zeros(length(pts_end),1);
%     for i=1:length(pts_start)
%         distances=dist_forward(i,:); % length: len(pts_end)
%         [~,idx]=sort(distances);
%         for j=1:length(distances)
%             if dist2Ant(pts_end(idx(j)))<dist2Ant(pts_start(i)) % partner is anterior neighbor
%                 startPt_forward(i)=pts_start(i);
%                 endPt_forward(i)=pts_end(idx(j));
%                 break
%             end
%         end
%     end
% 
%     startPt_backward=zeros(length(pts_end),1);
%     endPt_backward=zeros(length(pts_start),1);
%     for i=1:length(pts_end)
%         distances=dist_backward(i,:); % length: len(pts_start)
%         [~,idx]=sort(distances);
%         for j=1:length(distances)
%             if dist2Ant(pts_start(idx(j)))>dist2Ant(pts_end(i)) % partner is posterior neighbor
%                 startPt_backward(i)=pts_end(i);
%                 endPt_backward(i)=pts_start(idx(j));
%                 break
%             end
%         end
%     end
% 
%     if ~isempty(startPt_forward)
%         for i=1:length(startPt_forward)
%             if startPt_forward(i)~=0 && endPt_forward(i)~=0
%                 idx=find(startPt_backward==endPt_forward(i));
%                 if endPt_backward(idx)~=0 && endPt_backward(idx)==startPt_forward(i)
%                     p1=[x(startPt_forward(i)), y(startPt_forward(i)), z(startPt_forward(i))];
%                     p2=[x(endPt_forward(i)), y(endPt_forward(i)), z(endPt_forward(i))];
%                     if distance_pts(p1, p2)<=5
%                         startPt=[startPt,startPt_forward(i)];
%                         endPt=[endPt,endPt_forward(i)];
%                     end
%                 end
%             end
%         end
%         pts_start=setdiff(pts_start,startPt);
%         pts_end=setdiff(pts_end,endPt);
% 
%         dist_forward=zeros(length(pts_start),length(pts_end));
%         dist_Eu=zeros(length(pts_start),length(pts_end));
%         for i=1:length(pts_start)
%             for j=1:length(pts_end)
%                 pt=[x(pts_start(i)),y(pts_start(i)),z(pts_start(i))];
%                 pt1=[x(pts_end(j)),y(pts_end(j)),z(pts_end(j))];
%                 pt2=antPole;
%                 pt3=postPole;
%                 dist_Eu(i,j)=distance_pts(pt, pt1);
%                 dist_forward(i,j)=dist_Eu(i,j) + distance_pt2plane(pt1,pt,pt2,pt3);
%             end
%         end
% 
%         dist_backward=zeros(length(pts_end),length(pts_start));
%         for i=1:length(pts_end)
%             for j=1:length(pts_start)
%                 pt=[x(pts_end(i)),y(pts_end(i)),z(pts_end(i))];
%                 pt1=[x(pts_start(j)),y(pts_start(j)),z(pts_start(j))];
%                 pt2=antPole;
%                 pt3=postPole;
%                 dist_backward(i,j)=dist_Eu(i,j) + distance_pt2plane(pt1,pt,pt2,pt3);
%             end
%         end
%         currLen=length(pts_end);
%     end
% end
end
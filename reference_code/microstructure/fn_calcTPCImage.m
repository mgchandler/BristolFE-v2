function [r, wr] = fn_calcTPCImage(image, r, numrun)
% [r, wr] = calcTPC(model, r, numrun)
% Input parameters:
%     image  - segmented image, array
%     r      - correlation distances to be measured for, array 
%     numrun - number of random tests for each r, integer
% Output parameters
%     r      - correlation distances (same as input), array
%     wr     - correlation probability
% Program for Weixin Wang of Bristol University by 
% Ming Huang of Imperial College London, 19 November 2023
% https://www.sciencedirect.com/science/article/pii/S1359645403005469
% explains two point correlation

% preliminaries
xmin = [1,1];
xmax = size(image);
wr   = zeros(size(r));

% solve
for i = 1 : length(r)
    % first set of random points
    x1 = zeros(numrun, 2);
    for j = 1 : 2
        x1(:,j) = unifrnd(xmin(j), xmax(j), numrun, 1);
    end
    % second set of random points having distance r(i) to the first set
    angles = unifrnd(-pi,pi, numrun, 1); % random angles, probability of generation is equal
    x2  = x1 + r(i)*[cos(angles), sin(angles)];
    % round the points off to integral points
    x1 = round(x1);
    x2 = round(x2);
    % remove the points that are out of the model domain
    mm = zeros(numrun,1);
    mn = x2(:,1) - xmin(1);
    mx = x2(:,1) - xmax(1);
    mm(mn<0) = mm(mn<0) + 1;
    mm(mx>0) = mm(mx>0) + 1;
    mn = x2(:,2) - xmin(2);
    mx = x2(:,2) - xmax(2);
    mm(mn<0) = mm(mn<0) + 1;
    mm(mx>0) = mm(mx>0) + 1;
    x1 = x1(mm==0,:);
    x2 = x2(mm==0,:);
    % grain numbers
    flag1 = image(sub2ind(size(image),x1(:,1),x1(:,2)));
    flag2 = image(sub2ind(size(image),x2(:,1),x2(:,2)));
    % correlation
    wr(i) = sum(flag1==flag2)/numrun;
    disp([num2str(i), ' of ', num2str(length(r)), ' lengths done.'])
end

end
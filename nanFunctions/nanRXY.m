function [rxx, ryy, rxy, ryx] = nanRXY(X, Y)
% nanRXY - Compute cross covariance while ignoring data points with Nan
% assignment.
%[rxy,rxx,ryy,ryx] = nanRXY(X, Y)

numDims=size(X,2);
% Compute covariance with nan values
nanCov=nancov([X Y],'pairwise');

rxx=nanCov(1:numDims,1:numDims);
ryy=nanCov(numDims+1:end,numDims+1:end);
rxy=nanCov(1:numDims,numDims+1:end);
ryx=nanCov(numDims+1:end,1:numDims);

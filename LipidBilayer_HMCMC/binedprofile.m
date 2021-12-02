function varargout = binedprofile(varargin)
% BINEDPROFILE Reduce number of layers in electon density profile.
%   BINEDP = BINEDPROFILE(EDP) dynamically vary number of layers
%       for purpose of reducing number of layers. Consecutive layers whose
%       dipsersions vary by an amount of less than 1e-10 are binned. The
%       topmost (reference) and bottommost (substrate) are kept.
%   BINEDP = BINEDPROFILE(EDP,MINDIFF) takes specified minimum threshold
%       for the binning. If MINDIFF is not given, 1e-10 is used.
%
%   [BINEDP,IROWS] = BINEDPROFILE(EDP) returns index of row numbers that 
%       are kept from EDP.
%
% Zhang Jiang
% $Revision: 1.0 $  $Date: 2008/09/26 $
% $Revision: 1.1 $  $Date: 2011/06/28 $ Include IROWS output argument.
% $Revision: 1.2 $  $Date: 2015/12/03 $ Use gradient to bin.

%
% See also EDPROFILE

if nargin > 2
    error('Invalid input arguement ...');
end
edp = varargin{1};
minDiff = 1e-10;
if nargin == 2
    minDiff = varargin{2};
end

% --- new method using gradient
diffDelta = gradient(edp(:,2));
ind = find(abs(diffDelta)>=minDiff);
binedp = edp(ind,1:3);
%if binedp(1,1)~=edp(1,1)
    binedp = [edp(1,:); binedp];     % keep surface
    binedp(1) = binedp(1) - eps(edp(1,1));
%end 
%if binedp(end,1)~=edp(end,1)
    binedp = [binedp;edp(end,:)];     % keep substrate
    binedp(end,1) = binedp(end,1) + eps(edp(end,1));
%end 
ind = [1;ind;size(edp,1)];
varargout{1} = binedp;
if nargout == 2
    varargout{2} = ind;
end


% % --- old method
% diffDelta = abs(diff(edp(:,2)));
% ind = find(diffDelta<minDiff);
% %ind(1) = [];                % remove the topmost in order to keep the reference material
% diffInd = diff(ind);
% indInd = find(diffInd~=1);
% indInd = [0;indInd;length(ind)];
% binedp = edp;
% for iii = (length(indInd)-1):-1:1
%     binedp(ind(indInd(iii)+1)+1,2:end) = mean(binedp((ind(indInd(iii)+1)+1):(ind(indInd(iii+1))+1),2:end));
%     binedp((ind(indInd(iii)+1)+2):(ind(indInd(iii+1))+1),:) = [];
% end
% if binedp(end,1)~=edp(end,1)
%     binedp(end,:) = edp(end,:);     % keep substrate
% end 
% varargout{1} = binedp;
% if nargout == 2
%     [~,irows] = intersect(edp(:,1),binedp(:,1));
%     varargout{2} = irows;    
% end
%     



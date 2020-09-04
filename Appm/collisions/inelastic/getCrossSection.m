function [cs, csDataOut] = getCrossSection(e)
% get cross section at energy e

persistent csData

if isempty(csData)
    filename = 'cs_ArgonIonization_Rapp1965.csv';
    fprintf('Read data from file: %s \n', filename)
    csData = readmatrix(filename);
end

% interpolate data
cs = zeros(size(e));
if ~isempty(e)
    cs = interp1(csData(:,1), csData(:,2), e, 'spline', 'extrap');
end


% clip data outside of domain
idx_lo = e <= csData(1,1);
idx_hi = e >= csData(end,1);
cs( idx_lo ) = 0; 
cs( idx_hi ) = csData(end,2); 

csDataOut = csData;


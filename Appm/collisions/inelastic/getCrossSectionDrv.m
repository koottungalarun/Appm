function [ds,dCSdataOut] = getCrossSectionDrv(e)
% derivative of cross section data
persistent dCS

if isempty(dCS)
    [~,csData] = getCrossSection([]);
    
    n = 10000;
    e0 = [0 logspace(-3,3,n)];
    
    % numerical approximation of derivative
    de = 1e-2;
    csPlus = getCrossSection(e0 + de);
    csMinus = getCrossSection(e0 - de);
    diffNumApprox = (csPlus - csMinus) / (2*de);
    
    % ensure positive values at 'low' energies
    idx = e0 < 50;
    diffNumApprox(idx) = max(diffNumApprox(idx), 0);

    dCS(:,1) = e0;
    dCS(:,2) = diffNumApprox;
end

% interpolate data
ds = zeros(size(e));
if ~isempty(e)
    ds = interp1(dCS(:,1), dCS(:,2), e, 'spline', 'extrap');
end

% ensure positive values at low energies; this is needed because spline
% interpolation tries to approximate data smoothly, although there is a
% jump at ionization energy (E_i = 15.76 eV for Argon)
idx = e < 50;
ds(idx) = max(ds(idx), 0);

% clip data outside of domain
idx_lo = e <= dCS(1,1);
idx_hi = e >= dCS(end,1);
cs( idx_lo ) = 0; 
cs( idx_hi ) = 0; 

dCSdataOut = dCS;

end
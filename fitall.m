
% Fitting function for bilinear model
% Author: Ryan R Mahutga - Barocas Research Group - University of
% Minnesota, Minneapolis, MN, USA

% Last Update: 05/14/2020

% INPUTS: lam - Nx1 vector of tissue stretches
%       PKS - Nx1 vector of PK1 stresses
%       choice of lam_lo_lims or PK1_lo_lims - 2x1 values of small-strain limits
%       for fit in either stretch or PK1 stress
%       choice of lam_hi_lims or PK1_hi_lims - 2x1 values of large-strain
%       limits for fit in either stretch of PK1 stress

% OUTPUTS: SS_mod - small-strain modulus
%          lock_mod - large-strain modulus
%           lam_cross - transition stretch
%           PK!_cross - transition PK1 stress
%           p_low - polynomial fit for small-strain
%           p_high - polynomial fit for large-strain modulus
%           error - 2x1 values of fit errors for small-strain and
%           large-strain moduli

function [SS_mod, lock_mod, lam_cross, PK1_cross,  p_low, p_high, error] = fitall(lam, PKS, lam_lo_lims, PK1_lo_lims, lam_hi_lims, PK1_hi_lims)
        
        % determine if small-strain limits are stretch or PK1 stress
        if ~isempty(lam_lo_lims)
            idx = find(lam > lam_lo_lims(1) & lam<lam_lo_lims(2));
        elseif ~isempty(PK1_lo_lims)
            idx = find(PKS > PK1_lo_lims(1) & PKS<PK1_lo_lims(2));
        else
            errordlg('No limits defined or no data in specified limits');
        end
        
        % fit small-strain modulus
        [p_low, error(1)] = polyfit(lam(idx(1):idx(end)),PKS(idx(1):idx(end)),1);
        
        SS_mod = p_low(1);
        
        % determine if small-strain limits are stretch or PK1 stress
        if ~isempty(lam_hi_lims)
            idx2 = find(lam > lam_hi_lims(1) & lam<lam_hi_lims(2));
        elseif ~isempty(PK1_hi_lims)
            idx2 = find(PKS > PK1_hi_lims(1) & PKS<PK1_hi_lims(2));
        else
            errordlg('No limits defined or no data in specified limits');
        end
        
        % fit large-strain modulus
        [p_high, error(2)] = polyfit(lam(idx2(1):idx2(end)),PKS(idx2(1):idx2(end)),1);
        lock_mod = p_high(1) ;
       
        % determine transition stretch and stress at transition
       lam_cross = (p_high(2)-p_low(2))/(p_low(1)-p_high(1));
       PK1_cross = p_low(1)*lam_cross + p_low(2);

end
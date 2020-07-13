
% Solution for uniaxial stretch and transverse stretch
% Author: Ryan R Mahutga - Barocas Research Group - University of
% Minnesota, Minneapolis, MN, USA

% Last Update: 05/14/2020

% INPUTS: Ri - tissue inner radius (undeformed)
%         a - pin radius
%         H - tissue thickness (undeformed)
%         x0 - starting pin-to-pin distance (center to center)
%         delta - Nx1 displacement of the pin

% OUTPUTS: lam - Nx1 uniaxial tissue stretch for a given displacement
%          lam_t - Nx1 transverse tissue stretch for a given displacment

function [ lam, lam_t ] = lam_transverse(Ri,a, H, x0, delta )

mult=0.50; % thickness used for determining lambda and lambda_transverse 
%(0 is inner surface, 0.5 is centerline, and 1 is outer surface)

lam_g = 1.0; % initial guess for transverse lambda

L0=2*pi*(Ri + mult*H); % initial length

for k=1:length(delta)
    
    if k>1
        lam_g=lam_t(k-1); % update transverse stretch guess
    end
    
    Lcg = (2*pi()*(a+mult*lam_g*H)+2*(x0+delta(k))); % tissue current length at guessed transverse stretch
    
    lam(k) = Lcg/L0; % tissue stretch at guessed transverse stretch
    
    if Lcg>=L0 % if the tissue is stretched
        %eq = @(lam_p) (pi()*H*lam_p^3 + 2*(pi()*a+x0+delta(k))*lam_p^2 - pi()*(2*Ri+H))^2; % volume conservation equation
        
        p = [pi()*H, 2*(pi()*a+x0+delta(k)), 0, -pi()*(2*Ri+H)]; % coefficients of the third order polynomial above 
        
        rts = roots(p); % roots of above equation
        
       lam_t(k) = rts(rts<1.1 & rts>0); % transverse stretch
        
        Lc = (2*pi()*(a+mult*lam_t(k)*H)+2*(x0+delta(k))); % current length
        lam(k) = Lc/L0; % tissue stretch
    else % if the tissue is compressed
        lam_t(k)=lam_g;
    end
    
end

end


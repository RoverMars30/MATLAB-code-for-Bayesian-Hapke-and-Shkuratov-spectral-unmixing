function reff = Hapke_ssa2reff(ssa,i,e,g)
    %   角度单位是度

    mu0 = cosd(i);
    mu = cosd(e);
    
    N = size(ssa,2);    
    reff = zeros(size(ssa));
    % tmp1 = mu0/(mu0+mu);
    tmp1 = 1 / (mu0 + mu);
    for i = 1:N
        w = ssa(:,i);
        % tmp2 = (1+Back_fun(g))*Phase_fun(g)+H_fun(mu0,w).*H_fun(mu,w)-1;
        % tmp2 = (1+0)*1+H_fun(mu0,w).*H_fun(mu,w)-1;
        tmp2 = (1 + 0) * Phase_fun(g) + H_fun(mu0,w).*H_fun(mu,w) - 1;

        reff(:,i) = w/4*tmp1.*tmp2;
        
    end   
end

function B_g = Back_fun(g)
%     phi = 0.8;
%     h = -3/8*log(1-phi);
    h = 0.06;
    B_g = 1/(1+1/h*tand(g/2));   
end

function P_g = Phase_fun(g)
    b = -0.4; c = 0.25;   % Li Shuai 2011 JGR
    % b = 0.12; c = 0.6;  % Vincendon 2013 PSS
    P1_g = cosd(g); P2_g = 3/2*cosd(g)^2-1/2;
    P_g = 1+b*P1_g+c*P2_g;
end

function H = H_fun(mu,w)   
    r = sqrt(1-w);
    r0 = (1-r)./(1+r);    
    tmp = r0+(1-1/2*r0-r0*mu)*log((1+mu)/mu);
    
    H = 1./(1-(1-r)*mu.*tmp);
end
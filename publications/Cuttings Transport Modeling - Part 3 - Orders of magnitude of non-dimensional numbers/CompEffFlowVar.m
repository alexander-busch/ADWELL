function [ w_o, w_p, h_f, A_a_f, A_p_f, d_h_a_eff, d_h_p_eff ] = CompEffFlowVar( d_i, d_o, d_h, h_b )
%CompEffFlowVar Compute effective flow variables
% Returns
% - w_o = Outer width of annular bed
% - w_p = width of equivalent pipe bed
% - h_f = height of fluid above pipe bed
% - A_a_f = Fluids area in annulus
% - A_p_f = Fluids area in pipe
% - d_h_a_eff = Effective anular hydraulic diameter
% - d_h_a_eff = Effective pipe hydraulic diameter


% h_b=h_b(2)


% Eccentricity
e = 0;

% if length(d_h)==1
%     d_o = d_o.*ones(1,2);
%     d_i = d_i.*ones(1,2);
%     d_h = d_h.*ones(1,2);
%     r_o = d_o.*ones(1,2);
%     r_i = d_i.*ones(1,2);
%     r_h = d_h.*ones(1,2);
%     e = zeros(1,2);
% end



% Radii
r_o = d_o/2;
r_i = d_i/2;
r_h = d_h/2;

% Angle and outer width of annular bed
alpha = 2*acos(1-2.*h_b./d_o);
w_o = d_o.*sin(alpha/2);

% Preallocate annular bed area and wetted perimeter 
A_b=zeros(length(d_h),1);
P_f=zeros(length(d_h),1);

% Loop annular geometries (d_o, d_i, d_h)
for ii = 1:length(d_h)
    
    % Bed below drill pipe
    if h_b(ii) <= (r_o(ii)-e-r_i(ii))
        
        % Annular solids bed area A_b = f(h_b), Cayeux et al. (2014) 
        A_b(ii) = acos((r_o(ii)-h_b(ii))/r_o(ii))*r_o(ii)^2-(r_o(ii)-h_b(ii))*sqrt(r_o(ii)^2-(r_o(ii)-h_b(ii))^2);
        
        % Annular wetted perimeter
        P_f(ii) = pi.*d_o(ii)-alpha(ii).*r_o(ii) + w_o(ii) + pi.*d_i(ii)-alpha(ii).*r_i(ii);
        
    % Bed at drill pipe
    elseif h_b(ii) <= (r_o(ii)-e+r_i(ii))
        
        % Inner width of annular bed
        w_i = d_i(ii).*sin(alpha(ii)/2);
        
        % Annular solids bed area A_b = f(h_b), Cayeux et al. (2014) 
        A_b(ii) = acos((r_o(ii)-h_b(ii))/r_o(ii))*r_o(ii)^2-(r_o(ii)-h_b(ii))*sqrt(r_o(ii)^2-(r_o(ii)-h_b(ii))^2) - (acos((r_o(ii)-h_b(ii)-e)/r_i(ii))*r_i(ii)^2-(r_o(ii)-h_b(ii)-e)*sqrt(r_i(ii)^2+(r_o(ii)-h_b(ii)-e)^2));
                
        % Annular wetted perimeter
        P_f(ii) = pi.*d_o(ii)-alpha(ii).*r_o(ii) + w_o(ii)-w_i + pi.*d_i(ii)-alpha(ii).*r_i(ii) ;
        
    % Bed above drill pipe
    else
        
        % Annular solids bed area A_b = f(h_b), Cayeux et al. (2014) 
        A_b(ii) = acos((r_o(ii)-h_b(ii))/r_o(ii))*r_o(ii)^2-(r_o(ii)-h_b(ii))*sqrt(r_o(ii)^2-(r_o(ii)-h_b(ii))^2)-pi*r_i(ii)^2;
        
        % Annular wetted perimeter
        P_f(ii) = pi.*d_o(ii)-alpha(ii).*r_o(ii) + w_o(ii);
        
    end
end

% Fluids area in annulus
A_a_f = pi/4.*(d_o.^2-d_i.^2) - A_b;

% Effective hydraulic diameter for annulus
d_h_a_eff = 4.*A_a_f./P_f;

% Solids bed fraction in annulus (and pipe)
alpha_b = A_b./(pi/4.*(d_o.^2-d_i.^2));

% Solids bed area in pipe
A_b = (pi/4.*(d_h.^2)).*alpha_b;

% Fluids bed area in pipe
A_p_f = (pi./4*(d_h.^2))-A_b;

% Iterate angle alpha for pipe
alpha = 1;
alpha_old = 0;
eps = 0.0001;
while abs(alpha-alpha_old)>eps
    alpha_old = alpha;
    alpha = 8*A_b./d_h.^2+sin(alpha);
    % alpha*180/pi
end

% Solids bed height in pipe
h_b_p = r_h.*(1-cos(alpha./2));
% alpha_p = 2*acos(1-2.*h_b_p./d_h);

% Fluids height in pipe
h_f = d_h-h_b_p;

% Solids bed width in pipe
w_p = d_h.*sin(alpha./2);

% Effective hydraulic diameter for pipe
d_h_p_eff = 4.*A_p_f./(pi.*d_h-alpha.*r_h + w_p);

end


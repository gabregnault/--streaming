function [u_x_euler,u_y_euler,u_x_stokes,u_y_stokes] = streaming_0m_fct(rho,mu,R0,f0,a_1,x,y,X1,Y1,theta,r,L)

% Paramètres d'initialisation
m = 1;
nu = mu/rho; % viscosité cinématique

omega_1 = 2*pi*f0; % pulsation du mode 0
omega_m = omega_1; % pulsation du mode m imposée égale à omega_0

s_1 = a_1; % amplitude modale du mode 0

% Preliminary constant
delta_1 = sqrt(2*nu/omega_1); 
k_1 = (1+i)/delta_1;

% New variables and constants in the complex space
X = k_1*r;
x_1 = k_1*R0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF THE ANGLE AND RADIAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hankel functions + derivatives at x_1
H0_m = sqrt(pi/(2*x_1)) * besselh(m+1/2,1,x_1); % order m = 1
H0_m1 = sqrt(pi/(2*x_1)) * besselh(m+3/2,1,x_1); % order m+1
H1 = m/x_1*H0_m - H0_m1; % first derivative at the order m
H2 = (2/x_1)*H0_m1 - (1 - m*(m-1)/x_1^2)*H0_m; % second derivative at the order m
H3 = (m-2)/x_1*(m*(m-1)./x_1^2 - 1)*H0_m + (1 - (m^2+m+6)/x_1^2)*H0_m1; % third derivative at the order m

% Hankel functions + derivatives. General formulation
h0_m = @(X,m) sqrt(pi./(2.*X)) .* besselh(m+1/2,1,X); % order m = 1
h0_m1 = @(X,m) sqrt(pi./(2.*X)) .* besselh(m+3/2,1,X); % order m+1
h1 = @(X,m) m./X.*h0_m(X,m) - h0_m1(X,m); % first derivative at the order m
h2 = @(X,m) (2./X).*h0_m1(X,m) - (1 - m*(m-1)./X.^2).*h0_m(X,m); % second derivative at the order m
h3 = @(X,m) (m-2)./X.*(m*(m-1)./X.^2 - 1).*h0_m(X,m) + (1 - (m^2+m+6)./X.^2).*h0_m1(X,m); % third derivative at the order m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANT FOR THE FIRST ORDER VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b_m = 2*i*R0*(m+2)*omega_1*s_1 / ((m+1)*(-x_1^2*H2 + (m^2+3*m+2)*H0_m)); % from eq.19, paper part I
% a1 = 1/6*(k_1*R0)^2*H2*b_m;

b_m = -(3*i*R0*omega_1*s_1)/(x_1^2*H2 - 6*H0_m);
a1 = i*R0*omega_1*s_1*x_1^2*H2/(2*(x_1^2*H2 - 6*H0_m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS FOR THE EULERIAN AND STOKES VELOCITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = @(X) 1./X.^4.*(-x_1.^4.*H2.*conj(X.*h1(X,m) - h0_m(X,m)) - 6.*X.^3.*h1(X,m).*conj(h0_m(X,m)));
%%
% Constant C30 and C40
Q3 = integral(@(X)G(X).*X,x_1,k_1*15*R0,'Waypoints',k_1*linspace(R0,15*R0,1000)); % Integration up to "infinity"
Q4 = integral(@(X)G(X)./X,x_1,k_1*15*R0,'Waypoints',k_1*linspace(R0,15*R0,1000)); % Integration up to "infinity"


C30 = 1/30*Q3;
C40 = -1/70*Q4;

% Constant C10 and C20
% A = -C30*x_1^(5) - C40*x_1^(7) + x_1^(3)*H2*conj(2*H0_m + x_1*H1) - ...
%     6*x_1^2*H1*conj(H0_m);
% 
% B = -6*C30*x_1^(5) -16*C40*x_1^(7) + 6*x_1*H0_m*(6*conj(x_1)*conj(H1) -...
%     x_1^3*conj(H3)) + x_1^3*H2*(60*conj(H0_m) - 6*x_1*conj(H1) + ...
%     2*abs(x_1)^2*conj(H2) - x_1^3*conj(H3));

A = -C30*x_1^(5) - C40*x_1^(7) - x_1^(3)*H2*conj(2*H0_m + x_1*H1) - ...
    6*x_1^2*H1*conj(H0_m);

B = -6*C30*x_1^(5) -16*C40*x_1^(7) + x_1^(3)*conj(x_1)*conj(H3)*(x_1^2*H2 - 6*H0_m) + ...
    2*x_1^(3)*conj(H2)*(-2*x_1^2*H2 - 3*x_1*H1 + 6*H0_m) - ...
    12*x_1^2*conj(H0_m)*(4*x_1*H2 + 3*H1);

C10 = (B-6*A)/10;
C20 = (16*A-B)/(10*x_1^2);


for ii = 1:length(X1)
    for jj = 1:length(Y1)
        if theta(ii,jj) >= 0
            S = 0;
            S_prime = 0;
            if r(ii,jj) > 0.99*R0 % no calculation needed inside the bubble
                space = k_1*linspace(R0,r(ii,jj),1000); % imposed path for the integration

                q1 = integral(@(X) X.^6.*G(X),k_1*R0,X(ii,jj),'Waypoints',space);
                C1 = C10 - 1/70*q1;

                q2 = integral(@(X) X.^4.*G(X),k_1*R0,X(ii,jj),'Waypoints',space);
                C2 = C20 + 1/30*q2;

                q3 = integral(@(X) X.*G(X),k_1*R0,X(ii,jj),'Waypoints',space);
                C3 = C30 - 1/30*q3;

                q4 = integral(@(X) G(X)./X,k_1*R0,X(ii,jj),'Waypoints',space);
                C4 = C40 + 1/70*q4;
                
                F(ii,jj) = C1/X(ii,jj)^(3) + C2/X(ii,jj) + C3*X(ii,jj)^2 + C4*X(ii,jj)^4;
                F_prime(ii,jj) = -3*C1/X(ii,jj)^4 -  C2/X(ii,jj)^2 + 2*C3*X(ii,jj) + 4*C4*X(ii,jj)^3;
                
                
                % Eulerian velocities (eq.4 and eq.5)
                Vr_11(ii,jj) = - abs(b_m).^2./(6*nu.*r(ii,jj)).*(1-3*cos(theta(ii,jj)).^2).*real(F(ii,jj));
                Vt_11(ii,jj) = - abs(b_m).^2./(6*nu.*r(ii,jj)).*cos(theta(ii,jj)).*sqrt(1-cos(theta(ii,jj)).^2)*real(F(ii,jj) + X(ii,jj).*F_prime(ii,jj));
                
                % stokes velocities (eq.A19 and eq.A20)             
                Vsr_11(ii,jj) = abs(b_m).^2/(6*omega_1*r(ii,jj)^3)*(1-3*cos(theta(ii,jj)).^2)*real( 6*i*X(ii,jj)*conj(h0_m(X(ii,jj),m))*h1(X(ii,jj),m) + i*x_1^4*H2/X(ii,jj)^2*conj(2*h0_m(X(ii,jj),m) + X(ii,jj)*h1(X(ii,jj),m)) );
                Vst_11(ii,jj) = abs(b_m).^2/(6*omega_1*r(ii,jj)^3)*cos(theta(ii,jj))*sqrt(1 - cos(theta(ii,jj))^2)*real( 6*i*X(ii,jj)^2*h0_m(X(ii,jj),m)*conj(h2(X(ii,jj),m)) - i*x_1^4*H2/X(ii,jj)^2*conj(6*h0_m(X(ii,jj),m) - X(ii,jj)^2*h2(X(ii,jj),m)) );
                                
            else 
                F(ii,jj) = 0;
                F_prime(ii,jj) = 0;
                Vsr_11(ii,jj) = 0;
                Vst_11(ii,jj) = 0;
                Vr_11(ii,jj) = 0;
                Vt_11(ii,jj) = 0;
            end
        else 
            F(ii,jj) = 0;
            F_prime(ii,jj) = 0;
            Vsr_11(ii,jj) = 0;
            Vst_11(ii,jj) = 0;
            Vr_11(ii,jj) = 0;
            Vt_11(ii,jj) = 0;
        end
    end
end


% Use of symmetry to complet the map

for ii = 1:length(Y1)
    Vr_11(end-ii+1,:) = Vr_11(ii,:);
    Vt_11(end-ii+1,:) = -Vt_11(ii,:);
    Vsr_11(end-ii+1,:) = Vsr_11(ii,:);
    Vst_11(end-ii+1,:) = -Vst_11(ii,:);
end

% Vr_11(26,75)
% Vt_11(26,75)
% Vsr_11(26,75)
% Vst_11(26,75)


% Norm of the lagrangian velocity
VL11 = sqrt((Vsr_11+Vr_11).^2 + (Vst_11+Vt_11).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM POLAR TO CARTESIAN COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_x_euler = (Vr_11).*cos(theta) - (Vt_11).*sin(theta); 
u_y_euler = (Vr_11).*sin(theta) + (Vt_11).*cos(theta);

u_x_stokes = (Vsr_11).*cos(theta) - (Vst_11).*sin(theta); 
u_y_stokes = (Vsr_11).*sin(theta) + (Vst_11).*cos(theta);
end
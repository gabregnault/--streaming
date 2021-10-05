function [u_x_euler,u_y_euler,u_x_stokes,u_y_stokes] = streaming_0m_fct(m,rho,mu,R0,f0,a_1,a_m,delta_phi,x,y,X1,Y1,theta,r,L)

% Paramètres d'initialisation

nu = mu/rho; % viscosité cinématique

omega_1 = 2*pi*f0; % pulsation du mode 0
omega_m = omega_1; % pulsation du mode m imposée égale à omega_0

s_1 = a_1; % amplitude modale du mode 0
s_m = a_m*exp(-i*delta_phi); % amplitude du mode m




% Preliminary constant
delta_m = sqrt(2*nu/omega_m);
k_m = (1+i)/delta_m;

delta_1 = sqrt(2*nu/omega_1);
k_1 = (1+i)/delta_1;

% New variables and constants in the complex space
X_m = k_m*r;
X = k_1*r;
x_m0 = k_m*R0;
x_bar = k_1*R0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITION OF THE ANGLE AND RADIAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hankel functions + derivatives
% case m = 1
h1_0 = @(X) sqrt(pi./(2*X)) .* besselh(1+1/2,1,X); % order m
h1_0_m1 = @(X) sqrt(pi./(2*X)) .* besselh(1+3/2,1,X); % order m+1
h1_1 = @(X) 1./X.*h1_0(X) - h1_0_m1(X); % first derivative
h1_2 = @(X) 2./X.*h1_0_m1(X) - (1 - 1*(1-1)./X.^2).*h1_0(X); % second derivative
h1_3 = @(X) (1-2)./X.*(1*(1-1)./X.^2 - 1).*h1_0(X) + (1 - (1^2+1+6)./X.^2).*h1_0_m1(X); % third derivative

% case m > 1
hm_0 = @(X) sqrt(pi./(2*X)) .* besselh(m+1/2,1,X); % order m
hm_0_m1 = @(X) sqrt(pi./(2*X)) .* besselh(m+3/2,1,X); % order m+1
hm_1 = @(X) m./X.*hm_0(X) - hm_0_m1(X); % first derivative
hm_2 = @(X) 2./X.*hm_0_m1(X) - (1 - m*(m-1)./X.^2)*hm_0(X); % second derivative
hm_3 = @(X) (m-2)./X*(m*(m-1)./X.^2 - 1).*hm_0(X) + (1 - (m^2+m+6)./X.^2)*hm_0_m1(X); % third derivative

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANT FOR THE FIRST ORDER VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_m = 2*i*R0*(m+2)*omega_m*s_m / ((m+1)*(-x_bar^2*hm_2(x_bar) + (m^2+3*m+2)*hm_0(x_bar))); % from eq.19, paper part I
b_1 = 2*i*R0*(1+2)*omega_1*s_1 / ((1+1)*(-x_bar^2*h1_2(x_bar) + (1^2+3*1+2)*h1_0(x_bar))); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS FOR THE EULERIAN AND STOKES VELOCITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G1 = @(X) -x_bar^4./(6*X.^4)*conj(h1_2(x_bar)).*(2*hm_0(X) - X.*hm_1(X)) + 1./X.*conj(h1_0(X)).*hm_1(X) - 1./conj(X).*conj(h1_1(X)).*hm_0(X);
G2 = @(X) -m*(m+1)*hm_0(X).*conj(x_bar^4./(12*X.^4)*h1_2(x_bar) + 1./X.*h1_1(X)) + (m+1)*x_bar^(m+1)./(4*(m+2)*X.^(m+3)).*(x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar)).*conj(h1_0(X) - X.*h1_1(X));


f1 = @(X) (m*G1(X) - G2(X)).*X.^(-m);
f2 = @(X) ((m+1)*G1(X) + G2(X)).*X.^(4-m);
f3 = @(X) (3*G1(X) + 2*G2(X)).*X.^(2-m);
f4 = @(X) (2*m*(m+1)*G1(X) + G2(X)).*X.^(2-m);


space1 = k_1*linspace(R0,15*R0,1000); % for intergration up to "infinity"
I1 = integral(f1,x_bar,15*x_bar,'Waypoints',space1);
I2 = integral(f2,x_bar,15*x_bar,'Waypoints',space1);
I3 = integral(f3,x_bar,15*x_bar,'Waypoints',space1);
I4 = integral(f4,x_bar,15*x_bar,'Waypoints',space1);

gama_12_bar = -1/(2*(2*m+1)*(2*m+3)*(2*m+5))*I1;
gama_14_bar = 1/(2*(2*m+1)*(2*m-1)*(2*m-3))*I2;
gama_17_bar = -1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I3;
gama_57_bar = -1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I4;


S1 = @(X) (m+1)*( conj(3*x_bar^4./X.^4*h1_2(x_bar) + 6./X.*h1_1(X) - 6./X.^2.*h1_0(X)) .* ...
    ( (x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) - ...
    m*x_bar./X.*hm_0(X)) - conj(x_bar^4./X.^4.*h1_2(x_bar) - 6./X.^2.*h1_0(X)).* ...
    ( (x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/2 *(x_bar./X).^(m+2) + ...
    m*x_bar.*hm_1(X) - m*x_bar./X.*hm_0(X)) );

S2 = @(X) (m+1)/2*conj(x_bar^4./X.^4*h1_2(x_bar) + 6./X.*h1_1(X) + 6./X.^2.*h1_0(X)) .* ...
    ( (x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) - ...
    m*x_bar./X.*hm_0(X)) - conj(x_bar^4./X.^4.*h1_2(x_bar) - 6./X.^2.*h1_0(X)).* ...
    ( (x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) + ...
    x_bar.*hm_1(X) + x_bar./X.*hm_0(X));

S3 = @(X)  conj(2*x_bar^4./X.^4.*h1_2(x_bar) + 6./X.*h1_1(X)).* ...
    ((x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) + ...
    x_bar.*hm_1(X) + x_bar./X.*hm_0(X)) - conj(x_bar^4./X.^4.*h1_2(x_bar) - 6./X.^2.*h1_0(X)).* ...
    ((x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/2 *(x_bar./X).^(m+2) - ...
    x_bar.*X.*hm_2(X) - x_bar.*hm_1(X) + x_bar./X.*hm_0(X));

S4 = @(X) (m+1)/2*( conj(3*x_bar^4./X.^4*h1_2(x_bar) - 6*h1_2(X) - 6./X.*h1_1(X) + 6./X.^2.*h1_0(X)) .* ...
    ( (x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) - ...
    m*x_bar./X.*hm_0(X)) + conj(x_bar^4./X.^4*h1_2(x_bar) + 6./X.*h1_1(X) + 6./X.^2.*h1_0(X)) .* ...
    ( (m+1)*(x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) *(x_bar./X).^(m+2) + ...
    m*x_bar.*hm_1(X)) );


% S3_prime_xbar = conj(6/x_bar^2*h1_1(x_bar) - 14/x_bar*h1_2(x_bar)) * ...
%     ( (x_bar^2*hm_2(x_bar) - (m^2+m-2)*hm_0(x_bar))/(2*(m+2)) - x_bar*hm_1(x_bar) - hm_0(x_bar)) + ...
%     conj(2*h1_2(x_bar) - 6/x_bar*h1_1(x_bar)) * ( m*(m+1)/(2*x_bar)*hm_0(x_bar) - hm_1(x_bar) - 3/2*x_bar*hm_2(x_bar)) + ...
%     conj(4/x_bar*h1_2(x_bar) - 6/x_bar^2*h1_1(x_bar) + 12/x_bar^3*h1_0(x_bar)) * ...
%     (3/2*x_bar^2*hm_2(x_bar) + x_bar*hm_1(x_bar) - m*(m+1)/2*hm_0(x_bar)) - ...
%     conj(h1_2(x_bar) + 6/x_bar^2*h1_0(x_bar)) * ...
%     (x_bar^2*hm_3(x_bar) - (m-2)/2*x_bar*hm_2(x_bar) - hm_1(x_bar) + ((m+2)*(m^2+m-2)+2)/(2*x_bar)*hm_0(x_bar));
% 
% S4_prime_xbar = (m+1)/2 * ( conj(6/x_bar^3*h1_0(x_bar) - 6/x_bar^2*h1_1(x_bar) - 3/x_bar*h1_2(x_bar) + 3* h1_3(x_bar)) * ...
%     (x_bar^2*hm_2(x_bar)/(m+2) + (m+1)*hm_0(x_bar)) + ...
%     conj(9*h1_2(x_bar) + 6/x_bar*h1_1(x_bar) - 6/x_bar^2*h1_0(x_bar)) * ...
%     ( (m-2)*(m+1)/(2*x_bar)*hm_0(x_bar) + m*hm_1(x_bar) - x_bar/2*hm_2(x_bar) ) + ...
%     conj(6/x_bar^3*h1_0(x_bar) - 5/x_bar*h1_2(x_bar)) * ...
%     ( (m+1)*(x_bar^2*hm_2(x_bar) - (m^2+m-2)*hm_0(x_bar))/(m+2) - m*x_bar*hm_1(x_bar)) + ...
%     conj(h1_2(x_bar) - 6/x_bar*h1_1(x_bar) - 6/x_bar^2*h1_0(x_bar)) * ...
%     ( (m+1)*(m^2+m-2)*hm_0(x_bar)/(2*x_bar) - (3*m+1)/2*x_bar*hm_2(x_bar)));

S3_prime_xbar = conj(h1_2(x_bar) - 6/x_bar^2*h1_0(x_bar)) * ...
    ( x_bar^2*hm_3(x_bar) + (m+4)/2*x_bar*hm_2(x_bar) - 3*hm_1(x_bar) + (m+1)*(m^2+4*m-2)/(2*x_bar)*hm_0(x_bar) ) ...
    - conj(h1_2(x_bar) + 3/x_bar*h1_1(x_bar))*(x_bar*hm_2(x_bar)/(m+2) + 2*hm_1(x_bar) + (m+1)/x_bar*hm_0(x_bar) );
    

S4_prime_xbar = (m+1)/2 * (3*conj(x_bar*h1_3(x_bar) + 3*h1_2(x_bar) - 2/x_bar*h1_1(x_bar) + 2/x_bar^2*h1_0(x_bar))* ...
    ( (m+1)/x_bar*hm_0(x_bar) - x_bar*hm_2(x_bar)/(m+2)) + ...
    3/2*conj( h1_2(x_bar) + 2/x_bar*h1_1(x_bar) - 2/x_bar^2*h1_0(x_bar)) * ...
    ( x_bar*hm_2(x_bar) + 2*m*hm_1(x_bar) + (m^2 - m - 2)/x_bar*hm_0(x_bar) ) + ...
    2*conj(h1_2(x_bar) - 6/x_bar^2*h1_0(x_bar)) * ( (m+1)/(2*(m+2)*x_bar)*(x_bar^2*hm_2(x_bar) + (m^2+m-2)*hm_0(x_bar)) + m*hm_1(x_bar) ) + ...
    (m-1)/2 * conj( h1_2(x_bar) + 6/x_bar*h1_1(x_bar) + 6/x_bar^2*h1_0(x_bar) )* (x_bar*hm_2(x_bar) - (m+1)*(m+2)/x_bar*hm_0(x_bar) )  );


gama_11_bar = -gama_12_bar*x_bar^(2*m+5) + x_bar^(m+2)/(12*(2*m+1)*(2*m+3))* ...
    ( m*((m^2-3)*S1(x_bar) + (m^3-m+4)*S2(x_bar))/((m-1)*(m+2)) + m*(S3(x_bar) - x_bar*S3_prime_xbar) -...
    S4(x_bar) + x_bar*S4_prime_xbar);

gama_13_bar = -gama_14_bar*x_bar^(2*m-3) - x_bar^(m-2)/(12*(2*m-1)*(2*m+1))* ...
    ( (m+1)*((m^2+2*m-2)*S1(x_bar) - (m^3+3*m^2+2*m-4)*S2(x_bar))/((m-1)*(m+2)) + ...
    (m+1)*(S3(x_bar) - x_bar*S3_prime_xbar) + S4(x_bar) - x_bar*S4_prime_xbar);

gama_15_bar = -gama_17_bar*x_bar^(2*m+1) - x_bar^m/(12*(2*m-1)*(2*m+3))* ...
    ( m*(m+1)*( 3*S1(x_bar) + 2*(m^2+m-5)*S2(x_bar)) / ((m-1)*(m+2)) - ...
    3*(S3(x_bar) - x_bar*S3_prime_xbar) - 2*(S4(x_bar) - x_bar*S4_prime_xbar) );

gama_55_bar = -gama_57_bar*x_bar^(2*m+1) + x_bar^m/(12*(2*m-1)*(2*m+3))* ...
    ( ((2*m^4 + 4*m^3 - 9*m^2 - 11*m + 6)*S1(x_bar) + m*(m+1)*(3*m^2+3*m+2)*S2(x_bar))/((m-1)*(m+2)) + ...
    2*m*(m+1)*(S3(x_bar) - x_bar*S3_prime_xbar) + S4(x_bar) - x_bar*S4_prime_xbar );




Fct1 = @(X) (G2(X) - m*G1(X)).*X.^(m+5);
Fct2 = @(X) (m*G1(X) - G2(X)).*X.^(-m);
Fct3 = @(X) ((m+1)*G1(X) + G2(X)).*X.^(m+1);
Fct4 = @(X) ((m+1)*G1(X) + G2(X)).*X.^(4-m);
Fct5 = @(X) (3*G1(X) + 2*G2(X)).*X.^(m+3);
Fct6 = @(X) (3*G1(X) + 2*G2(X)).*X.^(2-m);
Fct7 = @(X) (2*m*(m+1)*G1(X) + G2(X)).*X.^(m+3);
Fct8 = @(X) (2*m*(m+1)*G1(X) + G2(X)).*X.^(2-m);


for ii = 1:length(X)
    for jj = 1:length(X)
        if theta(ii,jj) >= 0
            if r(ii,jj) > 0.999*R0
                space = k_1*linspace(R0,r(ii,jj),100); % integration line

                I_11 = integral(Fct1,x_bar,X(ii,jj),'Waypoints',space);
                gama_11 = gama_11_bar + 1/(2*(2*m+1)*(2*m+3)*(2*m+5))*I_11;

                I_12 = integral(Fct2,x_bar,X(ii,jj),'Waypoints',space);
                gama_12 = gama_12_bar + 1/(2*(2*m+1)*(2*m+3)*(2*m+5))*I_12;

                I_13 = integral(Fct3,x_bar,X(ii,jj),'Waypoints',space);
                gama_13 = gama_13_bar + 1/(2*(2*m+1)*(2*m-1)*(2*m-3))*I_13;

                I_14 = integral(Fct4,x_bar,X(ii,jj),'Waypoints',space);
                gama_14 = gama_14_bar - 1/(2*(2*m+1)*(2*m-1)*(2*m-3))*I_14;

                I_15 = integral(Fct5,x_bar,X(ii,jj),'Waypoints',space);
                gama_15 = gama_15_bar - 1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I_15;

                I_17 = integral(Fct6,x_bar,X(ii,jj),'Waypoints',space);
                gama_17 = gama_17_bar + 1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I_17;

                I_55 = integral(Fct7,x_bar,X(ii,jj),'Waypoints',space);
                gama_55 = gama_55_bar - 1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I_55;

                I_57 = integral(Fct8,x_bar,X(ii,jj),'Waypoints',space);
                gama_57 = gama_57_bar + 1/(2*(2*m-1)*(2*m+1)*(2*m+3))*I_57;

                F1 = gama_11*X(ii,jj)^(-m-2) + gama_12*X(ii,jj)^(m+3) + gama_13*X(ii,jj)^(2-m) + gama_14*X(ii,jj)^(m-1) + gama_15*X(ii,jj)^(-m) + gama_17*X(ii,jj)^(m+1);
                F2 = -(m+1)*gama_11*X(ii,jj)^(-m-2) - (m+1)*gama_12*X(ii,jj)^(m+3) + m*gama_13*X(ii,jj)^(2-m) + m*gama_14*X(ii,jj)^(m-1) + gama_55*X(ii,jj)^(-m) + gama_57*X(ii,jj)^(m+1);
                F1_prime = -(m+2)*gama_11*X(ii,jj)^(-m-3) + (m+3)*gama_12*X(ii,jj)^(m+2) + (2-m)*gama_13*X(ii,jj)^(1-m) + (m-1)*gama_14*X(ii,jj)^(m-2) - m*gama_15*X(ii,jj)^(-m-1) + (m+1)*gama_17*X(ii,jj)^m;
                F2_prime = (m+1)*(m+2)*gama_11*X(ii,jj)^(-m-3) - (m+1)*(m+3)*gama_12*X(ii,jj)^(m+2) + m*(2-m)*gama_13*X(ii,jj)^(1-m) + m*(m-1)*gama_14*X(ii,jj)^(m-2) - m*gama_55*X(ii,jj)^(-m-1) + (m+1)*gama_57*X(ii,jj)^m;

                % Legender functions
                P = legendre(m,cos(theta(ii,jj)));
                P1_m(ii,jj) = P(2);
                Pm(ii,jj) = P(1);

                % Eulerian velocities (eq.15 and eq.16)
                Vr_1m(ii,jj) = -1/nu/r(ii,jj)*real(conj(b_1)*b_m*( cos(theta(ii,jj))* Pm(ii,jj)*(m*(m+1)*F1 - 2*F2) + sqrt(1 - cos(theta(ii,jj))^2)*P1_m(ii,jj)*(F1 - F2)));
                Vt_1m(ii,jj) = -1/nu/r(ii,jj)*real(conj(b_1)*b_m*( cos(theta(ii,jj)) * P1_m(ii,jj)*(F1 + X(ii,jj)*F1_prime) + sqrt(1 - cos(theta(ii,jj))^2)*Pm(ii,jj)*(F2 + X(ii,jj)*F2_prime)));

                % Stokes velocities (eq.B52 and eq.B53)
                Vsr_1m(ii,jj) = -1/6/nu/R0*real(conj(b_1)*b_m*(cos(theta(ii,jj)).*Pm(ii,jj).*S1(X(ii,jj)) + sqrt(1 - cos(theta(ii,jj))^2).*P1_m(ii,jj).*S2(X(ii,jj))));
                Vst_1m(ii,jj) = -1/6/nu/R0*real(conj(b_1)*b_m*(cos(theta(ii,jj)).*P1_m(ii,jj).*S3(X(ii,jj)) + sqrt(1 - cos(theta(ii,jj))^2).*Pm(ii,jj).*S4(X(ii,jj))));
            else 
                Vr_1m(ii,jj) = 0;
                Vt_1m(ii,jj) = 0;
                Vsr_1m(ii,jj) = 0;
                Vst_1m(ii,jj) = 0;
            end
        else 
            Vr_1m(ii,jj) = 0;
            Vt_1m(ii,jj) = 0;
            Vsr_1m(ii,jj) = 0;
            Vst_1m(ii,jj) = 0;    
        end
    end
end

% Use of symmetry to complete the map
for ii = 1:length(Y1)
    Vr_1m(end-ii+1,:) = Vr_1m(ii,:);
    Vt_1m(end-ii+1,:) = -Vt_1m(ii,:);
    Vsr_1m(end-ii+1,:) = Vsr_1m(ii,:);
    Vst_1m(end-ii+1,:) = -Vst_1m(ii,:);
end

% Vr_1m(26,75)
% Vt_1m(26,75)
% Vsr_1m(26,75)
% Vst_1m(26,75)

% Norm of the lagrangian velocity
VL1m = sqrt((Vsr_1m+Vr_1m).^2 + (Vst_1m+Vt_1m).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM POLAR TO CARTESIAN COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% u_x = (Vsr_1m+Vr_1m).*cos(theta) - (Vst_1m+Vt_1m).*sin(theta); 
% u_y = (Vsr_1m+Vr_1m).*sin(theta) + (Vst_1m+Vt_1m).*cos(theta);


u_x_euler = (Vr_1m).*cos(theta) - (Vt_1m).*sin(theta); 
u_y_euler = (Vr_1m).*sin(theta) + (Vt_1m).*cos(theta);

u_x_stokes = (Vsr_1m).*cos(theta) - (Vst_1m).*sin(theta); 
u_y_stokes = (Vsr_1m).*sin(theta) + (Vst_1m).*cos(theta);

end
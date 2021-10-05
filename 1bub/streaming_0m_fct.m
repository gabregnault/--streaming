function [u_x_euler,u_y_euler,u_x_stokes,u_y_stokes] = streaming_0m_fct_V2(m,rho,mu,R0,f0,a_0,a_m,delta_phi,x,y,X1,Y1,theta,r,L)
% calcul de l'interaction du mode0 avec les autres modes oscillant à la
% même fréquence pour l'interface de calcul du streaming
format long % on prend en compte un max de décimale

% Paramètres d'initialisation

nu = 1e-6; % viscosité cinématique

omega_0 = 2*pi*f0; % pulsation du mode 0
omega_m = omega_0; % pulsation du mode m imposée égale à omega_0

% pour que theta varie entre 0 et pi et non -pi et pi

s_0 = a_0; % amplitude modale du mode 0
s_m = a_m*exp(-i*delta_phi); % amplitude du mode m


% calcul des constantes préliminaires
delta_m = sqrt(2*nu/omega_m); 
k_m = (1+i)/delta_m;

delta_0 = sqrt(2*nu/omega_0);
k_0 = (1+i)/delta_0;

% variables utilisées pour les calculs
X_m = k_m*r;
X = k_0*r;
x_m0 = k_m*R0;
x_00 = k_0*R0;

% calculs

% fonctions en hankel et ses dérivées en x_m0 (ou x_00)

% Fonctions de Hankel pour le mode 
h_0 = @(X,N) sqrt(pi./(2.*X)) .* besselh(N+1/2,1,X);
h_0_n1 = @(X,N) sqrt(pi./(2.*X)) .* besselh(N+3/2,1,X);
h_1 = @(X,N) N./X.*h_0(X,N) - h_0_n1(X,N); % dérivée première hankel m
h_2 = @(X,N) 2./X.*h_0_n1(X,N) - (1 - N*(N-1)./X.^2).*h_0(X,N);
h_3 = @(X,N) (N-2)./X.*(N*(N-1)./X.^2 - 1).*h_0(X,N) + (1 - (N^2+N+6)./X.^2).*h_0_n1(X,N);


a0 = i*R0*omega_0*s_0;
am = i*R0*omega_m*s_m*(-x_00^2*h_2(x_00,m) - (m-1)*(m+2)*h_0(x_00,m)) / ((m+1)*(-x_m0^2*h_2(x_00,m) + (m^2+3*m+2)*h_0(x_00,m)));
b_m = 2*i*R0*(m+2)*omega_m*s_m / ((m+1)*(-x_00^2*h_2(x_00,m) + (m^2+3*m+2)*h_0(x_00,m)));

am_bis = (-x_00^2*h_2(x_00,m) - (m-1)*(m+2)*h_0(x_00,m))/(2*(m+2))*b_m;

% definiton des fct nécessaires au calcul (pour les intégrales (87)->(90),
% (100) et (101).


G = @(X,M) (h_0(X,M) - X.*h_1(X,M))./X.^3;

step = 5*R0/7500;
alpha_inf = 15;
space_inf_x = k_m*(R0:step:alpha_inf*R0);

Int3 = trapz(space_inf_x,G(space_inf_x,m).*space_inf_x.^(3-m));
Int4 = trapz(space_inf_x,G(space_inf_x,m).*space_inf_x.^(1-m));

C30 = 1/(2*(2*m-1)*(2*m+1))*Int3;
C40 = -1/(2*(2*m+1)*(2*m+3))*Int4;

A_m = -C30*x_00^(2*m+1) - C40*x_00^(2*m+3) - x_00^(m-2)/2*((m+1)*h_0(x_00,m) + ...
    + 2*x_00*h_1(x_00,m) + x_00^2*h_2(x_00,m)/(m+2));

B_m = (1 - m^2)*C30*x_00^(2*m+1) - m*(m+2)*C40*x_00^(2*m+3) - x_00^(m-2)/2*( ...
    x_00^3*h_3(x_00,m) + (m^2 + 2*m +3)/(m+2)*x_00^2*h_2(x_00,m) + (m-1)*(m+2)*x_00*h_1(x_00,m) + ...
    + (m+1)*(m^2+4*m+1)*h_0(x_00,m));

C10 = (B_m - (m^2-1)*A_m)/(2*m+1);
C20 = (m*(m+2)*A_m - B_m)/((2*m+1)*x_00^2);


% calcul des termes de vitesse
Pm = zeros(L,L);
P1_m = zeros(L,L);
C_1 = zeros(L,L);
C_2 = zeros(L,L);
C_3 = zeros(L,L);
C_4 = zeros(L,L);

for ii = 1:length(X1)
    for jj = 1:length(Y1)
        if theta(ii,jj) >= 0
            S = 0;
            S_prime = 0;
            if r(ii,jj) > R0 % on ne calcul pas les vitesses à l'intérieur de la bulle ...
                space_x = k_m*(R0:step:r(ii,jj));
                
                I1 = trapz(space_x,G(space_x,m).*space_x.^(m+4));
                I2 = trapz(space_x,G(space_x,m).*space_x.^(m+2));
                I3 = trapz(space_x,G(space_x,m).*space_x.^(3-m));
                I4 = trapz(space_x,G(space_x,m).*space_x.^(1-m));
                
                C_1(ii,jj) = C10 - 1/(2*(2*m+1)*(2*m+3))*I1;
                C_2(ii,jj) = C20 + 1/(2*(2*m-1)*(2*m+1))*I2;
                C_3(ii,jj) = C30 - 1/(2*(2*m-1)*(2*m+1))*I3;
                C_4(ii,jj) = C40 + 1/(2*(2*m+1)*(2*m+3))*I4;
                               
                % associated legendre polynomial order 1
                P = legendre(m,cos(abs(theta(ii,jj))));
                P1_m(ii,jj) = P(2);
                Pm(ii,jj) = P(1);
                Vsr_0m(ii,jj) = m*(m+1)*R0*Pm(ii,jj)/2/omega_0/r(ii,jj)^4*real(i*conj(a0)*b_m*( (-x_00^2 * h_2(x_00,m) - (m-1)*(m+2)*h_0(x_00,m))/(2*(m+2))*(x_00/X(ii,jj))^(m+1) - h_0(X(ii,jj),m) - X(ii,jj)*h_1(X(ii,jj),m) )); % vitesse normale de stokes
                Vst_0m(ii,jj) = R0*P1_m(ii,jj)/2/omega_0/r(ii,jj)^4*real(i*conj(a0)*b_m*( 2*h_0(X(ii,jj),m) - X(ii,jj)^2*h_2(X(ii,jj),m) - (m+3)/(2*(m+2))*(-x_00^2*h_2(x_00,m) - (m-1)*(m+2)*h_0(x_00,m))*(x_00/X(ii,jj))^(m+1) ));
            else 
                P = legendre(m,cos(theta(ii,jj)));
                P1_m(ii,jj) = P(2);
                Pm(ii,jj) = P(1);
                Vsr_0m(ii,jj) = 0;
                Vst_0m(ii,jj) = 0;
            end
        else 
            P = legendre(m,cos(theta(ii,jj)));
            P1_m(ii,jj) = P(2);
            Vsr_0m(ii,jj) = 0;
            Vst_0m(ii,jj) = 0;
        end
    end
end

F = C_1./X.^(m+1) + C_2./X.^(m-1) + C_3.*X.^(m) + C_4.*X.^(m+2);
F_prime = -(m+1)*C_1./X.^(m+2) - (m-1)*C_2./X.^(m) + m*C_3.*X.^(m-1) + (m+2)*C_4.*X.^(m+1);

% calcul final pour les vitesse eulériennes
Vr_0m = m*(m+1)*R0./(2*nu.*r).*Pm.*real(k_0*conj(a0)*b_m.*F);  % vitesse normale Eulérienne
Vt_0m = R0./(2*nu.*r).*P1_m.*real(k_0*conj(a0)*b_m.*(F+X.*F_prime)); % vitesse tangentielle Eulérienne


% On complète l'espace avec la symétrie

for ii = 1:length(Y1)
    Vr_0m(end-ii+1,:) = Vr_0m(ii,:);
    Vt_0m(end-ii+1,:) = -Vt_0m(ii,:);
    Vsr_0m(end-ii+1,:) = Vsr_0m(ii,:);
    Vst_0m(end-ii+1,:) = -Vst_0m(ii,:);
end

% Vr_0m(26,75)
% Vt_0m(26,75)
% Vsr_0m(26,75)
% Vst_0m(26,75)

% passage en cartésien
u_x_euler = (Vr_0m).*cos(theta) - (Vt_0m).*sin(theta); % champ lointain
u_y_euler = (Vr_0m).*sin(theta) + (Vt_0m).*cos(theta);

u_x_stokes = (Vsr_0m).*cos(theta) - (Vst_0m).*sin(theta); % champ lointain
u_y_stokes = (Vsr_0m).*sin(theta) + (Vst_0m).*cos(theta);

end


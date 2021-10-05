function [u_x_euler,u_y_euler,u_x_stokes,u_y_stokes] = streaming_nm_fct(n,m,rho,mu,R0,f0,a_n,a_m,delta_phi,x,y,X1,Y1,theta,r,L)
% calcul de l'interaction de deux modes d'ordres n et m. n = m peut être
% considéré avec cette function.

nu = mu/rho; % viscosité cinématique

omega_n = 2*pi*f0;
omega_m = omega_n; % pulsation du mode m supposée égale à omega_m

s_n = a_n; % amplitude modale du mode 0
s_m = a_m*exp(-i*delta_phi); % amplitude du mode m


delta_m = sqrt(2*nu/omega_m);
k_m = (1+i)/delta_m;

delta_n = sqrt(2*nu/omega_n);
k_n = (1+i)/delta_n;

x_bar = k_n*R0;
x = k_m*r;

% calculs

% fonctions en hankel et ses dérivées en X
% Fonctions de Hankel pour le mode n
h_0 = @(X,N) sqrt(pi./(2.*X)) .* besselh(N+1/2,1,X);
h_0_n1 = @(X,N) sqrt(pi./(2.*X)) .* besselh(N+3/2,1,X);
h_1 = @(X,N) N./X.*h_0(X,N) - h_0_n1(X,N); % dérivée première hankel m
h_2 = @(X,N) 2./X.*h_0_n1(X,N) - (1 - N*(N-1)./X.^2).*h_0(X,N);
h_3 = @(X,N) (N-2)./X.*(N*(N-1)./X.^2 - 1).*h_0(X,N) + (1 - (N^2+N+6)./X.^2).*h_0_n1(X,N);



% équation (32) de la partie sur les mode 0 et m


b_m = 2*i*R0*(m+2)*omega_m*s_m / ((m+1)*(-x_bar^2*h_2(x_bar,m) + (m^2+3*m+2)*h_0(x_bar,m)));
a_m = i*R0*omega_m*s_m*(-x_bar^2*h_2(x_bar,m) - (m-1)*(m+2)*h_0(x_bar,m)) / ((m+1)*(-x_bar^2*h_2(x_bar,m) + (m^2+3*m+2)*h_0(x_bar,m)));

b_n = 2*i*R0*(n+2)*omega_n*s_n / ((n+1)*(-x_bar^2*h_2(x_bar,n) + (n^2+3*n+2)*h_0(x_bar,n)));
a_n = i*R0*omega_n*s_n*(-x_bar^2*h_2(x_bar,n) - (n-1)*(n+2)*h_0(x_bar,n)) / ((n+1)*(-x_bar^2*h_2(x_bar,n) + (n^2+3*n+2)*h_0(x_bar,n)));

% b_m = 2*i*R0*(m+2)*omega_m*s_m / ((m+1)*(x_bar^2*h_2(x_bar,m) + (m^2+3*m+2)*h_0(x_bar,m))); % constante équation (32)
% b_n = 2*i*R0*(n+2)*omega_n*s_n / ((n+1)*(x_bar^2*h_2(x_bar,n) + (n^2+3*n+2)*h_0(x_bar,n)));
% 
% a_n = i*R0*omega_n*s_n/(n+1) - n*h_0(x_bar,n)*b_n % eq15 papier 0-m
% a_m = i*R0*omega_m*s_m/(m+1) - m*h_0(x_bar,m)*b_m

% a_n = a_n;
% a_m = a_m;
% 
% b_n = b_n;
% b_m = b_m;


% Symboles de Wigner3j
alpha = @(k,n,m) - (2*k+1)/(k*(k+1))*sqrt((m+1)*m*(k+1)*k)*Wigner3j([n m k],[0 0 0])*Wigner3j([n m k],[0 1 -1]);
beta = @(k,n,m) (2*k+1)/(k*(k+1))*sqrt((n+1)*n*(k+1)*k*factorial(m+2)/factorial(m-2))*Wigner3j([n k m],[0 0 0])*Wigner3j([n k m],[1 1 -2]);
gamma = @(k,n,m) n*(n+1)/2*alpha(k,n,m) + m*(m+1)/2*alpha(k,m,n) - 1/2*beta(k,n,m) - 1/2*beta(k,m,n);

% fonctions G pour calcul (équation (51) et (52))
G1 = @(x,r) (n+1)./(2*nu*r.^2).*(k_n^2.*a_n.*conj(b_m).*(x_bar./x).^(n+1).*conj((n+1).*h_0(x,m) - x.*h_1(x,m)) - n.*k_n.^2.*b_n*conj(b_m).*(x.*h_1(x,n).*conj(h_0(x,m)) + conj(x).*h_0(x,n).*conj(h_1(x,m))));
G2 = @(x,r) (m+1)./(2*nu*r.^2).*(k_n^2*a_m*conj(b_n).*(x_bar./x).^(m+1).*conj((m+1)*h_0(x,n) - x.*h_1(x,n)) - m*k_n^2*b_m*conj(b_n)*(x.*h_1(x,m).*conj(h_0(x,n)) + conj(x).*h_0(x,m).*conj(h_1(x,n))));
G3 = @(x,r) -1./(2*nu*r.^2).*(conj(b_m)*k_n^2*conj(h_0(x,m)).*(a_n*(x_bar./x).^(n+1) - b_n*(h_0(x,n) + x.*h_1(x,n))) + conj(b_n)*k_n^2*conj(h_0(x,n)).*(a_m*(x_bar./x).^(m+1) - b_m*(h_0(x,m) + x.*h_1(x,m))));

Hk = @(k,n,m,x,r) alpha(k,n,m).*G1(x,r) + alpha(k,m,n).*G2(x,r) + gamma(k,n,m).*G3(x,r);

% Fonctions S1i 
S11 = @(x) -(m+1)./(x.^(2)).*a_n*conj(a_m).*(x_bar./x).^(n+m+3) + m*(m+1)*(x_bar./(x.^2))*b_n*conj(b_m).*h_1(x,n).*conj(h_0(x,m)) + a_m*conj(b_n)*(x_bar./x).^(m+2).*(m+1).*conj((m+1)./(x.^2).*h_0(x,n) + 1./x.*h_1(x,n));
S12 = @(x) -(n+1)./(x.^(2))*a_m*conj(a_n).*(x_bar./x).^(n+m+3) + n*(n+1)*(x_bar./(x.^2))*b_m*conj(b_n).*h_1(x,m).*conj(h_0(x,n)) + a_n*conj(b_m)*(x_bar./x).^(n+2).*(n+1).*conj((n+1)./(x.^2).*h_0(x,m) + 1./x.*h_1(x,m));
S13 = @(x) x_bar./(x.^(3))*b_n*conj(b_m).*h_0(x,n).*conj(h_0(x,m)) .* (n*(n+1) - m*(m+1)) + a_n*conj(b_m)*(x_bar./x).^(n+2).*conj(h_0(x,m))./(x.^2) * (m*(m+1) - n*(n+1)) + a_m*conj(b_n)*(x_bar./x).^(m+2).*conj(h_0(x,n))./(x.^2) * (n*(n+1) - m*(m+1));


% Fonctions S2i
S21 = @(x) x_bar./(x.^(3)) .* (a_n*(x_bar./x).^(n+1) + n*b_n*h_0(x,n)) .* conj( -(m+3)*a_m*(x_bar./x).^(m+1) + b_m*(2*h_0(x,m) - x.^2.*h_2(x,m)));
S22 = @(x) x_bar./(x.^(3)) .* (a_m*(x_bar./x).^(m+1) + m*b_m*h_0(x,m)) .* conj( -(n+3)*a_n*(x_bar./x).^(n+1) + b_n*(2*h_0(x,n) - x.^2.*h_2(x,n)));
S23 = @(x) x_bar./(x.^(3)) .* (a_n*(x_bar./x).^(n+1) - b_n*(h_0(x,n) + x.*h_1(x,n))) .* conj(a_m*(x_bar./x).^(m+1) - b_m*(h_0(x,m) + x.*h_1(x,m)));

% Dérivées des fonctions S2i
% S21_prime
A = @(x) x_bar./x.^3;
A_prime = @(x) -3*x_bar./x.^4;
B = @(x) a_n*(x_bar./x).^(n+1) + n*b_n*h_0(x,n);
B_prime = @(x) -a_n*(n+1)*x_bar^(n+1)./x.^(n+2) + n*b_n*h_1(x,n);
C = @(x) -(m+3)*a_m*(x_bar./x).^(m+1) + b_m*(2*h_0(x,m) - x.^2*h_2(x,m));
C_prime = @(x) (m+3)*(m+1)*a_m*x_bar^(m+1)./x.^(m+2) + b_m*(2*h_1(x,m) - 2*x.*h_2(x,m) - x.^2.*h_3(x,m));

S21_prime = A_prime(x_bar)*B(x_bar)*conj(C(x_bar)) + B_prime(x_bar)*A(x_bar)*conj(C(x_bar)) + conj(x_bar)./x_bar.^3*B(x_bar)*conj(C_prime(x_bar));

% S22_prime
B = @(x) a_m*(x_bar./x).^(m+1) + m*b_m*h_0(x,m);
B_prime = @(x) -a_m*(m+1)*x_bar^(m+1)./x.^(m+2) + m*b_m*h_1(x,m);
C = @(x) -(n+3)*a_n*(x_bar./x).^(n+1) + b_n*(2*h_0(x,n) - x.^2*h_2(x,n));
C_prime = @(x) (n+3)*(n+1)*a_n*x_bar^(n+1)./x.^(n+2) + b_n*(2*h_1(x,n) - 2*x.*h_2(x,n) - x.^2.*h_3(x,n));

S22_prime = A_prime(x_bar)*B(x_bar)*conj(C(x_bar)) + B_prime(x_bar)*A(x_bar)*conj(C(x_bar)) + conj(x_bar)./x_bar.^3*B(x_bar)*conj(C_prime(x_bar));

% S23_prime
D = @(x) a_n*(x_bar./x).^(n+1) - b_n*(h_0(x,n) + x.*h_1(x,n));
D_prime = @(x) -(n+1)*a_n*x_bar^(n+1)./x.^(n+2) - b_n*(2*h_1(x,n) + x.*h_1(x,n));
E = @(x) a_m*(x_bar./x).^(m+1) - b_m*(h_0(x,m) + x.*h_1(x,m));
E_prime = @(x) -(m+1)*a_m*x_bar^(m+1)./x.^(m+2) - b_m*(2*h_1(x,m) + x.*h_1(x,m));

S23_prime  = A_prime(x_bar)*D(x_bar)*conj(E(x_bar)) + D_prime(x_bar)*A(x_bar)*conj(E(x_bar)) + conj(x_bar)./x_bar.^3*D(x_bar)*conj(E_prime(x_bar));


% paramètres de calculs des intégrales
step = 5*R0/10000;
alpha_inf = 15;

% matrice de stockage des vitesses
Vr_nm = zeros(L,L);
Vt_nm = zeros(L,L);
Vsr_nm_bis = zeros(L,L);
Vst_nm_bis = zeros(L,L);

l_depart = n-m;
l_fin = n+m;
for k = (n-m:1:n+m)
    kindex = k+1;
%     fprintf('Boucle %d sur %d \n',k-l_depart+1,l_fin-l_depart+1);
    
    % définitions des paramètres nécessaires
    a_knm(kindex) = - (2*k+1)*sqrt((n+1)*n*(m+1)*m)*Wigner3j([k m n],[0 0 0])*Wigner3j([k m n],[0 1 -1]);
    b_knm(kindex) = (2*k+1)*Wigner3j([k n m],[0 0 0])*Wigner3j([k n m],[0 0 0]);

    % Fonctions en x_bar
    Tk_0(kindex) = - 1/(2*nu*R0)*(a_knm(kindex) - n*(n+1)*b_knm(kindex))*S11(x_bar) - 1/(2*nu*R0)*(a_knm(kindex) - m*(m+1)*b_knm(kindex))*S12(x_bar) + 1/(2*nu*R0)*a_knm(kindex)*S13(x_bar);
    Fk_0(kindex) = R0/(k*(k+1))*Tk_0(kindex);
    if k == 0
        Uk_0(kindex) = 0;
        Uk_0_prime(kindex) = 0;
    else
        Uk_0(kindex) = (n+1)/(2*nu*R0)*alpha(k,n,m)*S21(x_bar) + (m+1)/(2*nu*R0)*alpha(k,m,n)*S22(x_bar) + 1/(2*nu*R0)*S23(x_bar)*(m*(m+1)/2*alpha(k,m,n) - n*(n+1)/2*alpha(k,n,m) - 1/2*beta(k,n,m) + 1/2*beta(k,m,n));
        Uk_0_prime(kindex) = (n+1)/(2*nu*R0)*alpha(k,n,m)*S21_prime + (m+1)/(2*nu*R0)*alpha(k,m,n)*S22_prime + 1/(2*nu*R0)*S23_prime*(m*(m+1)/2*alpha(k,m,n) - n*(n+1)/2*alpha(k,n,m) - 1/2*beta(k,n,m) + 1/2*beta(k,m,n));
    end
    % début des calculs
    % Etape 1, calcul des Ck_i0
    space_inf_x = k_n*(R0:step:alpha_inf*R0);
    space_inf_r = (R0:step:alpha_inf*R0);
    
    if k == 0 
        Ck10(kindex) = 0;
        Ck20(kindex) = 0;
        Ck30(kindex) = 0;
        Ck40(kindex) = 0;
    else
        Ck10(kindex) = -1/(2*(2*k+1)*(2*k+3))*trapz(space_inf_x,space_inf_x.^(1-k)/(k_n^4).*Hk(k,n,m,space_inf_x,space_inf_r));
        Ck20(kindex) = 1/(2*(2*k-1)*(2*k+1))*trapz(space_inf_x,space_inf_x.^(3-k)/k_n^4.*Hk(k,n,m,space_inf_x,space_inf_r));
        Ck30(kindex) = - Ck20(kindex)*x_bar^(2*k-1) + R0/(2*(2*k+1))*x_bar^(k-1)*((k+3)/(k+1)*Tk_0(kindex) + Uk_0(kindex) - x_bar*Uk_0_prime(kindex));
        Ck40(kindex) = - Ck10(kindex)*x_bar^(2*k+3) - R0/(2*(2*k+1))*x_bar^(k+1)*((k-2)/k*Tk_0(kindex) + Uk_0(kindex) - x_bar*Uk_0_prime(kindex));
    end

    for ii = 1:length(x)
        for jj = 1:length(x)
            if theta(ii,jj) >= 0
                if r(ii,jj) > R0
                    % calcul des Cki
                    space_x  = k_n*(R0:step:r(ii,jj));
                    space_r  = (R0:step:r(ii,jj));
                    
                    % fonctions angulaires
                    Poly = legendre(k,cos(theta(ii,jj)));
                    data(kindex).Pk(ii,jj) = Poly(1);
                    
                    if k == 0
                        data(kindex).Pk_1(ii,jj) = 0; 
                        data(kindex).Ck1(ii,jj) = 0;
                        data(kindex).Ck2(ii,jj) = 0;
                        data(kindex).Ck3(ii,jj) = 0;
                        data(kindex).Ck4(ii,jj) = 0;
                    else
                        data(kindex).Pk_1(ii,jj) = Poly(2);
                        data(kindex).Ck1(ii,jj) = Ck10(kindex) + 1/(2*(2*k+1)*(2*k+3))*trapz(space_x,space_x.^(1-k)/k_n^4.*Hk(k,n,m,space_x,space_r));
                        data(kindex).Ck2(ii,jj) = Ck20(kindex) - 1/(2*(2*k-1)*(2*k+1))*trapz(space_x,space_x.^(3-k)/k_n^4.*Hk(k,n,m,space_x,space_r));
                        data(kindex).Ck3(ii,jj) = Ck30(kindex) + 1/(2*(2*k-1)*(2*k+1))*trapz(space_x,space_x.^(2+k)/k_n^4.*Hk(k,n,m,space_x,space_r));
                        data(kindex).Ck4(ii,jj) = Ck40(kindex) - 1/(2*(2*k+1)*(2*k+3))*trapz(space_x,space_x.^(4+k)/k_n^4.*Hk(k,n,m,space_x,space_r));
                    end                    
                end
            else 
                data(kindex).Ck1(ii,jj) = 0;
                data(kindex).Ck2(ii,jj) = 0;
                data(kindex).Ck3(ii,jj) = 0;
                data(kindex).Ck4(ii,jj) = 0;
                data(kindex).Pk(ii,jj) = 0;
                data(kindex).Pk_1(ii,jj) = 0;
            end
        end
    end
    Tk = - 1/(2*nu*R0)*(a_knm(kindex) - n*(n+1)*b_knm(kindex))*S11(x) - 1/(2*nu*R0)*(a_knm(kindex) - m*(m+1)*b_knm(kindex))*S12(x) + 1/(2*nu*R0)*a_knm(kindex)*S13(x);
    Uk = (n+1)/(2*nu*R0)*alpha(k,n,m)*S21(x) + (m+1)/(2*nu*R0)*alpha(k,m,n)*S22(x) + 1/(2*nu*R0)*S23(x)*(m*(m+1)/2*alpha(k,m,n) - n*(n+1)/2*alpha(k,n,m) - 1/2*beta(k,n,m) + 1/2*beta(k,m,n));

    
    data(kindex).Fk = data(kindex).Ck1.*x.^(k+2) + data(kindex).Ck2.*x.^(k) + data(kindex).Ck3./x.^(k-1) + data(kindex).Ck4./x.^(k+1);
    data(kindex).Fk_prime = (k+2)*data(kindex).Ck1.*x.^(k+1) + k*data(kindex).Ck2.*x.^(k-1) - (k-1)*data(kindex).Ck3.*x.^(-k) - (k+1)*data(kindex).Ck4.*x.^(-k-2);
    
    Vr_nm = Vr_nm + k*(k+1).*data(kindex).Fk.*data(kindex).Pk;
    Vt_nm = Vt_nm + (data(kindex).Fk + x.*data(kindex).Fk_prime).*data(kindex).Pk_1;
    Vsr_nm_bis = Vsr_nm_bis + real(Tk.*data(kindex).Pk);
    Vst_nm_bis = Vst_nm_bis + real(Uk.*data(kindex).Pk_1);
end

% calcul des Pn, Pm, Pn_1, Pm_1, Pn_1_prime et Pm_1_prime pour les vitesses
% de stokes

% initialisation
Pn = zeros(L,L);
Pm  = zeros(L,L);
Pm_1 = zeros(L,L);
Pn_1 = zeros(L,L);
Pn_1_prime = zeros(L,L);
Pm_1_prime = zeros(L,L);

for ii = 1:length(x)
    for jj = 1:length(x)
        if theta(ii,jj) >= 0
            if r(ii,jj) > R0
            % n
            Poly_n = legendre(n,cos(theta(ii,jj)));
            Poly_n_moins1 = legendre(n-1,cos(theta(ii,jj)));
            Pn(ii,jj) = Poly_n(1);
            Pn_1(ii,jj) = Poly_n(2);
            Pn_1_prime(ii,jj) = 1./(cos(theta(ii,jj))^2 - 1) * (cos(theta(ii,jj)) * n * Pn_1(ii,jj) - (n+1) * Poly_n_moins1(2));

            % m
            Poly_m = legendre(m,cos(theta(ii,jj)));
            Poly_m_moins1 = legendre(m-1,cos(theta(ii,jj)));
            Pm(ii,jj) = Poly_n(1);
            Pm_1(ii,jj) = Poly_m(2);
            Pm_1_prime(ii,jj) = 1./(cos(theta(ii,jj))^2 - 1) * (cos(theta(ii,jj)) * m * Pm_1(ii,jj) - (m+1) * Poly_m_moins1(2));
            end
        end
    end
end

Vr_nm = -1./r .* real(Vr_nm);
Vt_nm = -1./r .* real(Vt_nm);

Vsr_nm = -1/(2*nu*R0) * (Pn_1.*Pm_1 - n*(n+1)*Pn.*Pm).*real(S11(x)) ...
    - 1/(2*nu*R0) * (Pn_1.*Pm_1 - m*(m+1)*Pn.*Pm).*real(S12(x)) ...
    + 1/(2*nu*R0)*Pn_1.*Pm_1.*real(S13(x));

Vst_nm  = (n+1)/(2*nu*R0)*Pn.*Pm_1.*real(S21(x)) ...
    + (m+1)/(2*nu*R0)*Pm.*Pn_1.*real(S22(x)) ...
    + 1/(2*nu*R0)*sqrt(1-cos(theta).^2).*(Pn_1.*Pm_1_prime - Pm_1.*Pn_1_prime).*real(S23(x));


%%
for ii = 1:length(Y1)
    Vr_nm(end-ii+1,:) = Vr_nm(ii,:);
    Vt_nm(end-ii+1,:) = -Vt_nm(ii,:);
    Vsr_nm(end-ii+1,:) = Vsr_nm(ii,:);
    Vst_nm(end-ii+1,:) = -Vst_nm(ii,:);
end


% passage en coordonnées cartésiennes

u_x_euler = (Vr_nm).*cos(theta) - (Vt_nm).*sin(theta); % champ lointain
u_y_euler = (Vr_nm).*sin(theta) + (Vt_nm).*cos(theta);

u_x_stokes = (Vsr_nm).*cos(theta) - (Vst_nm).*sin(theta); % champ lointain
u_y_stokes = (Vsr_nm).*sin(theta) + (Vst_nm).*cos(theta);


end


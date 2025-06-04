%%  Application de la SVD : compression d'images

clear all
close all

% Lecture de l'image
I = imread('BD_Asterix_Colored.jpg');
I = rgb2gray(I);
I = double(I);

[q, p] = size(I)

% Décomposition par SVD
fprintf('Décomposition en valeurs singulières\ain')
tic
[U, S, V] = svd(I);
toc

l = min(p,q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On choisit de ne considérer que 200 vecteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 200 vecteurs utilisés pour la reconstruction et on affiche l'image tous les 40 vecteurs (pas)
inter = 1:40:(200+40);
inter(end) = 200;

% vecteur pour stocker la différence entre l'image et l'image reconstruite

differenceSVD = zeros(size(inter,2), 1);

% images reconstruites en utilisant de 1 à 200 vecteurs
ti = 0;
td = 0;
for k = inter

    % Calcul de l'image de rang k
    Im_k = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';

    % Affichage de l'image reconstruite
    ti = ti+1;
    figure(ti)
    colormap('gray')
    imagesc(Im_k), axis equal
    
    % Calcul de la différence entre les 2 images RMSE 
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    pause
end

% Figure des différences entre l'image réelle et les images reconstruites
ti = ti+1;
figure(ti)
hold on 
plot(inter, differenceSVD, 'rx')
ylabel('RMSE')
xlabel('rank k')
pause


% Plugger les différentes méthodes : eig, puissance itérée et les 4 versions de la "subspace iteration method" 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUELQUES VALEURS PAR DÉFAUT DE PARAMÈTRES, 
% VALEURS QUE VOUS DEVEZ FAIRE ÉVOLUER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolérance
eps = 1e-8;
% nombre d'itérations max pour atteindre la convergence
maxit = 10000;

% taille de l'espace de recherche (m)
search_space = 400;

% pourcentage que l'on se fixe
percentage = 0.99;

% p pour les versions 2 et 3 (attention p déjà utilisé comme taille)
puiss = 1;

%%%%%%%%%%%%%
% À COMPLÉTER
%%%%%%%%%%%%%

%%
% calcul des couples propres et calcul des valeurs singulières
%%

%%% avec la methode eig
M = I*I';
tic
[U_eig, D_eig] = eig(M);
t_eig = toc;
[d, indice] = sort(diag(D_eig), 'descend');
U_eig = U_eig(:, indice);
sigma_eig = sqrt(d);

% avec la methode power_v11
tic
[Vp11, Dp11, n_ev_p11, itv_p11, flag_p11] = power_v11(M, search_space, percentage, eps, maxit);
t_p11 = toc;
[d_p11, ip11] = sort(diag(Dp11), 'descend');
Vp11 = Vp11(:, ip11);
sigma_p11 = sqrt(d_p11);

% avec la methode power_v12
tic
[Vp12, Dp12, n_ev_p12, itv_p12, flag_p12] = power_v12(M, search_space, percentage, eps, maxit);
t_p12 = toc;
[d_p12, ip12] = sort(diag(Dp12), 'descend');
Vp12 = Vp12(:, ip12);
sigma_p12= sqrt(d_p12);

% avec la methode subspace_iter_v0
tic
[U0, D0, it0, flag0] = subspace_iter_v0(M, search_space, eps, maxit);
t_sub0 = toc;
[d0, i0] = sort(diag(D0), 'descend');
U0 = U0(:, i0);
sigma_subv0 = sqrt(d0);

% avec la methode subspace_iter_v1
tic
[U1, D1, it1, flag1] = subspace_iter_v1(M, search_space, percentage, eps, maxit);
t_sub1 = toc;
[d1, i1] = sort(diag(D1), 'descend');
U1 = U1(:, i1);
sigma_subv1 = sqrt(d1);

% avec la methode subspace_iter_v2
tic
[U2, D2, it2, flag2] = subspace_iter_v2(M, search_space, percentage, puiss, eps, maxit);
t_sub2 = toc;
[d2, i2] = sort(diag(D2), 'descend');
U2 = U2(:, i2);
sigma_subv2 = sqrt(d2);

% avec la methode subspace_iter_v3
tic
[U3, D3, it3, flag3] = subspace_iter_v3(M, search_space, percentage, puiss, eps, maxit);
t_sub3 = toc;
[d3, i3] = sort(diag(D3), 'descend');
U3 = U3(:, i3);
sigma_subv3 = sqrt(d3);
%%
% calcul de l'autre ensemble de vecteurs
%%
V_eig = (I' * U_eig ) ./ sigma_eig';
V11 = (I' * Vp11) ./ sigma_p11';
V12 = (I' * Vp12) ./ sigma_p12';
V_sub0 = (I' * U0) ./ sigma_subv0';
V_sub1 = (I' * U1) ./ sigma_subv1';
V_sub2 = (I' * U2) ./ sigma_subv2';
V_sub3 = (I' * U3) ./ sigma_subv3';

% TODO
%%
% calcul des meilleures approximations de rang faible
%%
rmse_eig  = zeros(length(inter),1);
rmse_v11  = zeros(length(inter),1);
rmse_v12  = zeros(length(inter),1);
rmse_subv0 = zeros(length(inter),1);
rmse_subv1 = zeros(length(inter),1);
rmse_subv2 = zeros(length(inter),1);
rmse_subv3 = zeros(length(inter),1);

for i = 1:length(inter)
    k = inter(i);

    Ieig = U_eig(:,1:k) * diag(sigma_eig(1:k)) * V_eig(:,1:k)';
    
    k_v11 = min(k, size(Vp11,2));
    powerv11 = Vp11(:,1:k_v11) * diag(sigma_p11(1:k_v11)) * V11(:,1:k_v11)';
    
    k_v12 = min(k, size(Vp12,2));
    powerv12 = Vp12(:,1:k_v12) * diag(sigma_p12(1:k_v12)) * V12(:,1:k_v12)';
    
    k_subv0 = min(k, size(U0,2));
    subspace0 = U0(:,1:k_subv0) * diag(sigma_subv0(1:k_subv0)) * V_sub0(:,1:k_subv0)';
    
    k_subv1 = min(k, size(U1,2));
    subspace1 = U1(:,1:k_subv1) * diag(sigma_subv1(1:k_subv1)) * V_sub1(:,1:k_subv1)';
    
    k_subv2 = min(k, size(U2,2));
    subspace2 = U2(:,1:k_subv2) * diag(sigma_subv2(1:k_subv2)) * V_sub2(:,1:k_subv2)';
    
    k_subv3 = min(k, size(U3,2));
    subspace3 = U3(:,1:k_subv3) * diag(sigma_subv3(1:k_subv3)) * V_sub3(:,1:k_subv3)';
    
    rmse_eig(i) = sqrt(mean((I(:)-Ieig(:)).^2));
    rmse_v11(i) = sqrt(mean((I(:)-powerv11(:)).^2));
    rmse_v12(i) = sqrt(mean((I(:)-powerv12(:)).^2));
    rmse_subv0(i) = sqrt(mean((I(:)-subspace0(:)).^2));
    rmse_subv1(i) = sqrt(mean((I(:)-subspace1(:)).^2));
    rmse_subv2(i) = sqrt(mean((I(:)-subspace2(:)).^2));
    rmse_subv3(i) = sqrt(mean((I(:)-subspace3(:)).^2));
end

for k = 1:40:(200+40)
    % power_v11
    k_bis = min(k, size(Vp11,2));
    Im_k = Vp11(:,1:k_bis) * diag(sigma_p11(1:k_bis)) * V11(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);
    imagesc(Im_k), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) power\\_v11 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % power_v12
    k_bis = min(k, size(Vp12,2));
    Im_k = Vp12(:,1:k_bis) * diag(sigma_p12(1:k_bis)) * V12(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);
    imagesc(Im_k), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) power\\_v12 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v0
    k_bis = min(k, size(U0,2));
    Im_k = U0(:,1:k_bis) * diag(sigma_subv0(1:k_bis)) * V_sub0(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);      
    imagesc(Im_k), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) subspace\\_iter\\_v0 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v1
    k_bis = min(k, size(U1,2));
    Im_k = U1(:,1:k_bis) * diag(sigma_subv1(1:k_bis)) * V_sub1(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);
    imagesc(Im_k), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) subspace\\_iter\\_v1 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v2
    k_bis = min(k, size(U2,2));
    Im_k = U2(:,1:k_bis) * diag(sigma_subv2(1:k_bis)) * V_sub2(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);
    imagesc(Im_k), axis equal
    td = td +  1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) subspace\\_iter\\_v2 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v3
    k_bis = min(k, size(U3,2));
    Im_k = U3(:,1:k_bis) * diag(sigma_subv3(1:k_bis)) * V_sub3(:,1:k_bis)';
    ti = ti+1;
    figure(ti)
    colormap('gray')
    Im_k = real(Im_k);
    imagesc(Im_k), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    fprintf('Différence (RMSE) subspace\\_iter\\_v3 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end


figure;
plot(inter, rmse_eig, '-o');
title('RMSE vs. Rank k - eig');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_v11, '-x');
title('RMSE vs. Rank k - power\_v11');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_v12, '-^');
title('RMSE vs. Rank k - power\_v12');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_subv0, '-s');
title('RMSE vs. Rank k - subspace\_v0');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_subv1, '-d');
title('RMSE vs. Rank k - subspace\_v1');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_subv2, '-v');
title('RMSE vs. Rank k - subspace\_v2');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;

figure;
plot(inter, rmse_subv3, '-p');
title('RMSE vs. Rank k - subspace\_v3');
xlabel('rank k'); 
ylabel('RMSE'); 
grid on;


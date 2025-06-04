%%  Application de la SVD : compression d'images couleur

clear all
close all

% Lecture de l'image
I = imread('BD_Asterix_Colored.jpg');
I = double(I);

[q, p, c] = size(I)

% Séparation des canaux RGB
Ir = I(:,:,1);
Ig = I(:,:,2);
Ib = I(:,:,3);

% Décomposition par SVD pour chaque canal
fprintf('Décomposition en valeurs singulières\n')
tic
[Ur, Sr, Vr] = svd(Ir);
[Ug, Sg, Vg] = svd(Ig);
[Ub, Sb, Vb] = svd(Ib);
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

    % Calcul de l'image de rang k pour chaque canal
    Im_rk = Ur(:, 1:k)*Sr(1:k, 1:k)*Vr(:, 1:k)';
    Im_gk = Ug(:, 1:k)*Sg(1:k, 1:k)*Vg(:, 1:k)';
    Im_bk = Ub(:, 1:k)*Sb(1:k, 1:k)*Vb(:, 1:k)';
    
    % Assemblage de l'image couleur
    Im_k = cat(3, Im_rk, Im_gk, Im_bk);

    % Affichage de l'image reconstruite
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    title(['Reconstruction avec k = ' num2str(k)])
    
    % Calcul de la différence entre les 2 images (RMSE : Root Mean Square Error)
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    pause
end

% Figure des différences entre l'image réelle et les images reconstruites
ti = ti+1;
figure(ti)
hold on 
plot(inter, differenceSVD, 'rx')
ylabel('RMSE')
xlabel('rank k')
title('Erreur de reconstruction (SVD)')
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
Mr = Ir*Ir';
Mv = Ig*Ig';
Iv = Ig;
Mb = Ib*Ib';
% pour chaque canal indépendamment

% Rouge
tic
[Ur_eig, Dr_eig] = eig(Mr);
t_eig_r = toc;
[dr, indice_r] = sort(diag(Dr_eig), 'descend');
Ur_eig = Ur_eig(:, indice_r);
sigma_r_eig = sqrt(dr);

% Vert
tic
[Ug_eig, Dg_eig] = eig(Mv);
t_eig_g = toc;
[dg, indice_g] = sort(diag(Dg_eig), 'descend');
Ug_eig = Ug_eig(:, indice_g);
sigma_g_eig = sqrt(dg);

% Bleu
tic
[Ub_eig, Db_eig] = eig(Mb);
t_eig_b = toc;
[db, indice_b] = sort(diag(Db_eig), 'descend');
Ub_eig = Ub_eig(:, indice_b);
sigma_b_eig = sqrt(db);

% avec la methode power_v11
% Rouge
tic
[Vp11r, Dp11r, n_ev_p11r, itv_p11r, flag_p11r] = power_v11(Mr, search_space, percentage, eps, maxit);
t_p11r = toc;
[d_p11r, ip11r] = sort(diag(Dp11r), 'descend');
Vp11r = Vp11r(:, ip11r);
sigma_p11r = sqrt(d_p11r);

% Vert
tic
[Vp11v, Dp11v, n_ev_p11v, itv_p11v, flag_p11v] = power_v11(Mv, search_space, percentage, eps, maxit);
t_p11v = toc;
[d_p11v, ip11v] = sort(diag(Dp11v), 'descend');
Vp11v = Vp11v(:, ip11v);
sigma_p11v = sqrt(d_p11v);

% Bleu
tic
[Vp11b, Dp11b, n_ev_p11b, itv_p11b, flag_p11b] = power_v11(Mb, search_space, percentage, eps, maxit);
t_p11b = toc;
[d_p11b, ip11b] = sort(diag(Dp11b), 'descend');
Vp11b = Vp11b(:, ip11b);
sigma_p11b = sqrt(d_p11b);


% avec la methode power_v12
% Rouge
tic
[Vp12r, Dp12r, n_ev_p12r, itv_p12r, flag_p12r] = power_v12(Mr, search_space, percentage, eps, maxit);
t_p12r = toc;
[d_p12r, ip12r] = sort(diag(Dp12r), 'descend');
Vp12r = Vp12r(:, ip12r);
sigma_p12r = sqrt(d_p12r);

% Vert
tic
[Vp12v, Dp12v, n_ev_p12v, itv_p12v, flag_p12v] = power_v12(Mv, search_space, percentage, eps, maxit);
t_p12v = toc;
[d_p12v, ip12v] = sort(diag(Dp12v), 'descend');
Vp12v = Vp12v(:, ip12v);
sigma_p12v = sqrt(d_p12v);

% Bleu
tic
[Vp12b, Dp12b, n_ev_p12b, itv_p12b, flag_p12b] = power_v12(Mb, search_space, percentage, eps, maxit);
t_p12b = toc;
[d_p12b, ip12b] = sort(diag(Dp12b), 'descend');
Vp12b = Vp12b(:, ip12b);
sigma_p12b = sqrt(d_p12b);


% avec la methode subspace_iter_v0
% Rouge
tic
[U0r, D0r, it0r, flag0r] = subspace_iter_v0(Mr, search_space, eps, maxit);
t_sub0r = toc;
[d0r, i0r] = sort(diag(D0r), 'descend');
U0r = U0r(:, i0r);
sigma_subv0r = sqrt(d0r);

% Vert
tic
[U0v, D0v, it0v, flag0v] = subspace_iter_v0(Mv, search_space, eps, maxit);
t_sub0v = toc;
[d0v, i0v] = sort(diag(D0v), 'descend');
U0v = U0v(:, i0v);
sigma_subv0v = sqrt(d0v);

% Bleu
tic
[U0b, D0b, it0b, flag0b] = subspace_iter_v0(Mb, search_space, eps, maxit);
t_sub0b = toc;
[d0b, i0b] = sort(diag(D0b), 'descend');
U0b = U0b(:, i0b);
sigma_subv0b = sqrt(d0b);


% avec la methode subspace_iter_v1
% Rouge
tic
[U1r, D1r, it1r, flag1r] = subspace_iter_v1(Mr, search_space, percentage, eps, maxit);
t_sub1r = toc;
[d1r, i1r] = sort(diag(D1r), 'descend');
U1r = U1r(:, i1r);
sigma_subv1r = sqrt(d1r);

% Vert
tic
[U1v, D1v, it1v, flag1v] = subspace_iter_v1(Mv, search_space, percentage, eps, maxit);
t_sub1v = toc;
[d1v, i1v] = sort(diag(D1v), 'descend');
U1v = U1v(:, i1v);
sigma_subv1v = sqrt(d1v);

% Bleu
tic
[U1b, D1b, it1b, flag1b] = subspace_iter_v1(Mb, search_space, percentage, eps, maxit);
t_sub1b = toc;
[d1b, i1b] = sort(diag(D1b), 'descend');
U1b = U1b(:, i1b);
sigma_subv1b = sqrt(d1b);

% avec la methode subspace_iter_v2
% Rouge
tic
[U2r, D2r, it2r, flag2r] = subspace_iter_v2(Mr, search_space, percentage, puiss, eps, maxit);
t_sub2r = toc;
[d2r, i2r] = sort(diag(D2r), 'descend');
U2r = U2r(:, i2r);
sigma_subv2r = sqrt(d2r);

% Vert
tic
[U2v, D2v, it2v, flag2v] = subspace_iter_v2(Mv, search_space, percentage, puiss, eps, maxit);
t_sub2v = toc;
[d2v, i2v] = sort(diag(D2v), 'descend');
U2v = U2v(:, i2v);
sigma_subv2v = sqrt(d2v);

% Bleu
tic
[U2b, D2b, it2b, flag2b] = subspace_iter_v2(Mb, search_space, percentage, puiss, eps, maxit);
t_sub2b = toc;
[d2b, i2b] = sort(diag(D2b), 'descend');
U2b = U2b(:, i2b);
sigma_subv2b = sqrt(d2b);


% avec la methode subspace_iter_v3
% Rouge
tic
[U3r, D3r, it3r, flag3r] = subspace_iter_v3(Mr, search_space, percentage, puiss, eps, maxit);
t_sub3r = toc;
[d3r, i3r] = sort(diag(D3r), 'descend');
U3r = U3r(:, i3r);
sigma_subv3r = sqrt(d3r);

% Vert
tic
[U3v, D3v, it3v, flag3v] = subspace_iter_v3(Mv, search_space, percentage, puiss, eps, maxit);
t_sub3v = toc;
[d3v, i3v] = sort(diag(D3v), 'descend');
U3v = U3v(:, i3v);
sigma_subv3v = sqrt(d3v);

% Bleu
tic
[U3b, D3b, it3b, flag3b] = subspace_iter_v3(Mb, search_space, percentage, puiss, eps, maxit);
t_sub3b = toc;
[d3b, i3b] = sort(diag(D3b), 'descend');
U3b = U3b(:, i3b);
sigma_subv3b = sqrt(d3b);


%%
% calcul de l'autre ensemble de vecteurs
%%
% V pour eig
Vr_eig = (Ir' * Ur_eig) ./ sigma_r_eig';
Vg_eig = (Ig' * Ug_eig) ./ sigma_g_eig';
Vb_eig = (Ib' * Ub_eig) ./ sigma_b_eig';

% V pour power_v11
V11r = (Ir' * Vp11r) ./ sigma_p11r';
V11v = (Iv' * Vp11v) ./ sigma_p11v';
V11b = (Ib' * Vp11b) ./ sigma_p11b';

% V pour power_v12
V12r = (Ir' * Vp12r) ./ sigma_p12r';
V12v = (Ig' * Vp12v) ./ sigma_p12v';
V12b = (Ib' * Vp12b) ./ sigma_p12b';

% V pour subspace_iter_v0
V0r = (Ir' * U0r) ./ sigma_subv0r';
V0v = (Ig' * U0v) ./ sigma_subv0v';
V0b = (Ib' * U0b) ./ sigma_subv0b';

% V pour subspace_iter_v1
V1r = (Ir' * U1r) ./ sigma_subv1r';
V1v = (Ig' * U1v) ./ sigma_subv1v';
V1b = (Ib' * U1b) ./ sigma_subv1b';

% V pour subspace_iter_v2
V2r = (Ir' * U2r) ./ sigma_subv2r';
V2v = (Ig' * U2v) ./ sigma_subv2v';
V2b = (Ib' * U2b) ./ sigma_subv2b';

% V pour subspace_iter_v3
V3r = (Ir' * U3r) ./ sigma_subv3r';
V3v = (Ig' * U3v) ./ sigma_subv3v';
V3b = (Ib' * U3b) ./ sigma_subv3b';

%%
% calcul des meilleures approximations de rang faible
%%

rmse_eig = zeros(max(inter),1);
rmse_v11 = zeros(max(inter),1);
rmse_v12  = zeros(max(inter),1);
rmse_subv0 = zeros(max(inter),1);
rmse_subv1 = zeros(max(inter),1);
rmse_subv2 = zeros(max(inter),1);
rmse_subv3 = zeros(max(inter),1);

for k = inter
    % eig 
    Irk = Ur_eig(:,1:k) * diag(sigma_r_eig(1:k)) * Vr_eig(:,1:k)';
    Igk = Ug_eig(:,1:k) * diag(sigma_g_eig(1:k)) * Vg_eig(:,1:k)';
    Ibk = Ub_eig(:,1:k) * diag(sigma_b_eig(1:k)) * Vb_eig(:,1:k)';
    Irec_k = cat(3, Irk, Igk, Ibk);
    rmse_eig(k) = sqrt(mean((I(:) - Irec_k(:)).^2));

    % power_v11
    k_v11r = min(k, size(Vp11r,2));
    powerv11r = Vp11r(:,1:k_v11r) * diag(sigma_p11r(1:k_v11r)) * V11r(:,1:k_v11r)';
    k_v11v = min(k, size(Vp11v,2));
    powerv11v = Vp11v(:,1:k_v11v) * diag(sigma_p11v(1:k_v11v)) * V11v(:,1:k_v11v)';
    k_v11b = min(k, size(Vp11b,2));
    powerv11b = Vp11b(:,1:k_v11b) * diag(sigma_p11b(1:k_v11b)) * V11b(:,1:k_v11b)';
    Iv11_k = cat(3, powerv11r, powerv11v, powerv11b);
    rmse_v11(k) = sqrt(mean((I(:)-Iv11_k(:)).^2));

    % power_v12
    k_v12r = min(k, size(Vp12r,2));
    powerv12r = Vp12r(:,1:k_v12r) * diag(sigma_p12r(1:k_v12r)) * V12r(:,1:k_v12r)';
    k_v12v = min(k, size(Vp12v,2));
    powerv12v = Vp12v(:,1:k_v12v) * diag(sigma_p12v(1:k_v12v)) * V12v(:,1:k_v12v)';
    k_v12b = min(k, size(Vp12b,2));
    powerv12b = Vp12b(:,1:k_v12b) * diag(sigma_p12b(1:k_v12b)) * V12b(:,1:k_v12b)';
    Iv12_k = cat(3, powerv12r, powerv12v, powerv12b);
    rmse_v12(k) = sqrt(mean((I(:)-Iv12_k(:)).^2));

    % subspace_iter_v0
    k_subv0r = min(k, size(U0r,2));
    subv0r = U0r(:,1:k_subv0r) * diag(sigma_subv0r(1:k_subv0r)) * V0r(:,1:k_subv0r)';
    k_subv0v = min(k, size(U0v,2));
    subv0v = U0v(:,1:k_subv0v) * diag(sigma_subv0v(1:k_subv0v)) * V0v(:,1:k_subv0v)';
    k_subv0b = min(k, size(U0b,2));
    subv0b = U0b(:,1:k_subv0b) * diag(sigma_subv0b(1:k_subv0b)) * V0b(:,1:k_subv0b)';
    Isubv0_k = cat(3, subv0r, subv0v, subv0b);
    rmse_subv0(k) = sqrt(mean((I(:)-Isubv0_k(:)).^2));

    % subspace_iter_v1
    k_subv1r = min(k, size(U1r,2));
    subv1r = U1r(:,1:k_subv1r) * diag(sigma_subv1r(1:k_subv1r)) * V1r(:,1:k_subv1r)';
    k_subv1v = min(k, size(U1v,2));
    subv1v = U1v(:,1:k_subv1v) * diag(sigma_subv1v(1:k_subv1v)) * V1v(:,1:k_subv1v)';
    k_subv1b = min(k, size(U1b,2));
    subv1b = U1b(:,1:k_subv1b) * diag(sigma_subv1b(1:k_subv1b)) * V1b(:,1:k_subv1b)';
    Isubv1_k = cat(3, subv1r, subv1v, subv1b);
    rmse_subv1(k) = sqrt(mean((I(:)-Isubv1_k(:)).^2));

    % subspace_iter_v2
    k_subv2r = min(k, size(U2r,2));
    subv2r = U2r(:,1:k_subv2r) * diag(sigma_subv2r(1:k_subv2r)) * V2r(:,1:k_subv2r)';
    k_subv2v = min(k, size(U2v,2));
    subv2v = U2v(:,1:k_subv2v) * diag(sigma_subv2v(1:k_subv2v)) * V2v(:,1:k_subv2v)';
    k_subv2b = min(k, size(U2b,2));
    subv2b = U2b(:,1:k_subv2b) * diag(sigma_subv2b(1:k_subv2b)) * V2b(:,1:k_subv2b)';
    Isubv2_k = cat(3, subv2r, subv2v, subv2b);
    rmse_subv2(k) = sqrt(mean((I(:)-Isubv2_k(:)).^2));

    % subspace_iter_v3
    k_subv3r = min(k, size(U3r,2));
    subv3r = U3r(:,1:k_subv3r) * diag(sigma_subv3r(1:k_subv3r)) * V3r(:,1:k_subv3r)';
    k_subv3v = min(k, size(U3v,2));
    subv3v = U3v(:,1:k_subv3v) * diag(sigma_subv3v(1:k_subv3v)) * V3v(:,1:k_subv3v)';
    k_subv3b = min(k, size(U3b,2));
    subv3b = U3b(:,1:k_subv3b) * diag(sigma_subv3b(1:k_subv3b)) * V3b(:,1:k_subv3b)';
    Isubv3_k = cat(3, subv3r, subv3v, subv3b);
    rmse_subv3(k) = sqrt(mean((I(:)-Isubv3_k(:)).^2));

end

figure;
plot(inter, rmse_eig(inter), '-o',    ...
     inter, rmse_v11(inter), '-x',    ...
     inter, rmse_v12(inter), '-^',    ...
     inter, rmse_subv0(inter),  '-s',    ...
     inter, rmse_subv1(inter),  '-d',    ...
     inter, rmse_subv2(inter),  '-v',    ...
     inter, rmse_subv3(inter),  '-p');
legend('eig','power\_v11','power\_v12','v0','v1','v2','v3','Location','northeast');
xlabel('rank k');
ylabel('RMSE');
title('Reconstruction RMSE vs. k for all methods');
grid on;



%% affichage couleur pour toutes les méthodes
for k = 1:40:(200+40)
    % eig
    k_r = min(k, size(Ur_eig,2));
    Irk = Ur_eig(:,1:k_r) * diag(sigma_r_eig(1:k_r)) * Vr_eig(:,1:k_r)';
    k_g = min(k, size(Ug_eig,2));
    Igk = Ug_eig(:,1:k_g) * diag(sigma_g_eig(1:k_g)) * Vg_eig(:,1:k_g)';
    k_b = min(k, size(Ub_eig,2));
    Ibk = Ub_eig(:,1:k_b) * diag(sigma_b_eig(1:k_b)) * Vb_eig(:,1:k_b)';
    Im_k = cat(3, Irk, Igk, Ibk);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) eig pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % power_v11
    k_r = min(k, size(Vp11r,2));
    powerv11r = Vp11r(:,1:k_r) * diag(sigma_p11r(1:k_r)) * V11r(:,1:k_r)';
    k_g = min(k, size(Vp11v,2));
    powerv11v = Vp11v(:,1:k_g) * diag(sigma_p11v(1:k_g)) * V11v(:,1:k_g)';
    k_b = min(k, size(Vp11b,2));
    powerv11b = Vp11b(:,1:k_b) * diag(sigma_p11b(1:k_b)) * V11b(:,1:k_b)';
    Im_k = cat(3, powerv11r, powerv11v, powerv11b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) power\\_v11 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % power_v12
    k_r = min(k, size(Vp12r,2));
    powerv12r = Vp12r(:,1:k_r) * diag(sigma_p12r(1:k_r)) * V12r(:,1:k_r)';
    k_g = min(k, size(Vp12v,2));
    powerv12v = Vp12v(:,1:k_g) * diag(sigma_p12v(1:k_g)) * V12v(:,1:k_g)';
    k_b = min(k, size(Vp12b,2));
    powerv12b = Vp12b(:,1:k_b) * diag(sigma_p12b(1:k_b)) * V12b(:,1:k_b)';
    Im_k = cat(3, powerv12r, powerv12v, powerv12b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) power\\_v12 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v0
    k_r = min(k, size(U0r,2));
    subv0r = U0r(:,1:k_r) * diag(sigma_subv0r(1:k_r)) * V0r(:,1:k_r)';
    k_g = min(k, size(U0v,2));
    subv0v = U0v(:,1:k_g) * diag(sigma_subv0v(1:k_g)) * V0v(:,1:k_g)';
    k_b = min(k, size(U0b,2));
    subv0b = U0b(:,1:k_b) * diag(sigma_subv0b(1:k_b)) * V0b(:,1:k_b)';
    Im_k = cat(3, subv0r, subv0v, subv0b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) subspace\\_iter\\_v0 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v1
    k_r = min(k, size(U1r,2));
    subv1r = U1r(:,1:k_r) * diag(sigma_subv1r(1:k_r)) * V1r(:,1:k_r)';
    k_g = min(k, size(U1v,2));
    subv1v = U1v(:,1:k_g) * diag(sigma_subv1v(1:k_g)) * V1v(:,1:k_g)';
    k_b = min(k, size(U1b,2));
    subv1b = U1b(:,1:k_b) * diag(sigma_subv1b(1:k_b)) * V1b(:,1:k_b)';
    Im_k = cat(3, subv1r, subv1v, subv1b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) subspace\\_iter\\_v1 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v2
    k_r = min(k, size(U2r,2));
    subv2r = U2r(:,1:k_r) * diag(sigma_subv2r(1:k_r)) * V2r(:,1:k_r)';
    k_g = min(k, size(U2v,2));
    subv2v = U2v(:,1:k_g) * diag(sigma_subv2v(1:k_g)) * V2v(:,1:k_g)';
    k_b = min(k, size(U2b,2));
    subv2b = U2b(:,1:k_b) * diag(sigma_subv2b(1:k_b)) * V2b(:,1:k_b)';
    Im_k = cat(3, subv2r, subv2v, subv2b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) subspace\\_iter\\_v2 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

for k = 1:40:(200+40)
    % subspace_iter_v3
    k_r = min(k, size(U3r,2));
    subv3r = U3r(:,1:k_r) * diag(sigma_subv3r(1:k_r)) * V3r(:,1:k_r)';
    k_g = min(k, size(U3v,2));
    subv3v = U3v(:,1:k_g) * diag(sigma_subv3v(1:k_g)) * V3v(:,1:k_g)';
    k_b = min(k, size(U3b,2));
    subv3b = U3b(:,1:k_b) * diag(sigma_subv3b(1:k_b)) * V3b(:,1:k_b)';
    Im_k = cat(3, subv3r, subv3v, subv3b);
    ti = ti+1;
    figure(ti)
    imagesc(uint8(Im_k)), axis equal
    td = td + 1;
    differenceSVD(td) = sqrt(mean((I(:)-Im_k(:)).^2));
    fprintf('Différence (RMSE) subspace\\_iter\\_v3 pour k = %d : %f\n', k, differenceSVD(td));
    pause
end

figure;
plot(inter, rmse_eig(inter), '-o',    ...
     inter, rmse_v11(inter), '-x',    ...
     inter, rmse_v12(inter), '-^',    ...
     inter, rmse_subv0(inter),  '-s',    ...
     inter, rmse_subv1(inter),  '-d',    ...
     inter, rmse_subv2(inter),  '-v',    ...
     inter, rmse_subv3(inter),  '-p');
legend('eig','power\_v11','power\_v12','v0','v1','v2','v3','Location','northeast');
xlabel('rank k');
ylabel('RMSE');
title('Reconstruction RMSE vs. k for all methods');
grid on;

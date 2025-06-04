load('A_100_1.mat'); A1 = A;
load('A_100_2.mat'); A2 = A;
load('A_100_3.mat'); A3 = A;
load('A_100_4.mat'); A4 = A;

D1 = eig(A1);
D2 = eig(A2);
D3 = eig(A3);
D4 = eig(A4);

figure;
subplot(2,2,1);
stem(sort(D1, 'descend'), 'filled'); title('Matrix A\_100\_1');
xlabel('Index'); ylabel('Valeur Propre');

subplot(2,2,2);
stem(sort(D2, 'descend'), 'filled'); title('Matrix A\_100\_2');
xlabel('Index'); ylabel('Valeur Propre');

subplot(2,2,3);
stem(sort(D3, 'descend'), 'filled'); title('Matrix A\_100\_3');
xlabel('Index'); ylabel('Valeur Propre');

subplot(2,2,4);
stem(sort(D4, 'descend'), 'filled'); title('Matrix A\_100\_4');
xlabel('Index'); ylabel('Valeur Propre');

sgtitle('Distribution des Valeur Propre pour les 4 types de matrices');

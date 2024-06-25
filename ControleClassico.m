clear 
clc
close all

% Definir os parâmetros do sistema
k_p = 105e3;  % N/m
k_m = 30.25e3;  % N/m
k_cinto = 52.5e3;  % N/m
K = 21.678;  % Nm/rad
k_b = 40e3;  % N/m
b_b = 220; % Ns/m
b_p = 325;  % Ns/m
b_m = 2408.04;  % Ns/m
b_cinto = 250;  % Ns/m
B_rot = 2.4692;  % Nms/rad
l = 0.19;  % m
M_total = 70;  % kg
M_t = 42.609;  % kg
M_m = 21.238;  % kg
M_c = 6.153;  % kg
J_c = 0.02;  % kgm^2
M_b = 35;  % kg
g = 9.81; % m/s^2

% Definindo as matrizes
n1 = ((1 + M_c * l^2) / (J_c + M_c * l^2) / (M_c + M_m))^(-1);
n2 = (-M_c * l / (J_c + M_c * l^2));
n3 = (M_c + M_m + ((M_c^2) * (l^2)) / (J_c + M_c * l^2));
A = [
    0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 1;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)
];

B = [
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    1 / M_b, k_b / M_b, b_b / M_b;
    0, 0, 0;
    0, k_p / M_m, b_p / M_m;
    0, 0, 0
];

C =  [
    1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0;
    -(k_b + k_cinto) / M_b, k_cinto / M_b, 0, 0, -(b_b + b_cinto) / M_b, b_cinto / M_b, 0, 0;
    k_cinto * n3, -(k_cinto + k_m) * n3, k_m * n3, M_c * l * (M_c * l * g - K) / (J_c + M_c * l^2) * n3, b_cinto / (M_t + M_c), -(b_cinto + b_m) / (M_t + M_c), b_m / (M_t + M_c), -B_rot * n3 / (J_c + M_c * l^2);
    0, k_m / M_m, -(k_m + k_p) / M_m, 0, 0, b_m / M_m, -(b_m + b_p) / M_m, 0;
    n1 * n2 * k_cinto, n1 * n2 * (-k_cinto - k_m), k_m * n1 * n2, (M_c * l * g - K) / (J_c + M_c * l^2), n1 * n2 * b_cinto, (-b_cinto - b_m) * n1 * n2, b_m * n1 * n2, -B_rot / (J_c + M_c * l^2)
];

D =  [
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    1 / M_b, k_b / M_b, b_b / M_b;
    0, 0, 0;
    0, k_p / M_m, b_p / M_m;
    0, 0, 0
];

% Definindo o espaço de estados
EspacoDeEstados = ss(A, B, C, D);
FT = tf(EspacoDeEstados(5, 1));
G = FT;

% Traçar o lugar das raízes variando Kp
figure(1)
rlocus(G);  % Corrigido para traçar o lugar das raízes da planta G
grid on
title('Lugar das Raízes para Variação de Kp');

% Ganho proporcional ótimo
Kp = 0.45;

% Função de transferência em malha fechada com Kp
FT_MF = feedback(Kp * G, 1);

% Função de transferência integrativa
s = tf('s');
FT_I = FT_MF / s;

figure(2)
rlocus(FT_I);  % Traçar lugar das raízes do sistema com ação integrativa
grid on 
title('Lugar das Raízes para Variação de Ki');

% Ganho integral ótimo
Ki = 20;

% Função de transferência em malha fechada com Kp e Ki
FT_PI = feedback(Kp * G + Ki / s * G, 1);  % Ajustado para FT_PI

% Função de transferência derivativa
FT_d = s * FT_PI;

figure(3)
rlocus(FT_d);  % Traçar lugar das raízes do sistema com ação derivativa
grid on 
title('Lugar das Raízes para Variação de Kd');

% Ganho derivativo
Kd = 0.0182;

% Controlador PID
PID = pid(Kp, Ki, 0);

% Função de transferência em malha fechada com o controlador PID
T = feedback(PID * G, 1);

figure(4)
margin(T);

tempoDeSimulacao = 2; %s
amostra = 0.0001;
tempo = 0:amostra:tempoDeSimulacao;
x0 = [0; 0; 0; 0; 0; 0; 0; 0];

u = 55438.083333333*ones(size(tempo));


[y, t] = lsim(G, u, tempo, x0);
[y_pi, t_pi] = lsim(T, u, tempo, x0);
y = y/9.81;
y_pi = y_pi/9.81;

figure(5)
plot(t, y, t_pi, y_pi)
xlabel('Tempo (s)');
ylabel('Aceleração (Gs)');
title('Aceleraçao do Banco');
legend('Malha aberta', 'PI')
% hold on;
% % Coordenadas do retângulo (área a ser ampliada)
% x1 = 0.05;
% x2 = 0.25;
% y1 = 25;
% y2 = 75;
% 
% % Desenhar o retângulo
% rectangle('Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', 'r', 'LineWidth', 2);
% 
% % Definir as coordenadas do novo eixo (axes) para o gráfico ampliado
% axes('Position', [0.6, 0.6, 0.25, 0.25]);
% box on;
% 
% % Desenhar o gráfico ampliado
% plot(t, y, t_pi, y_pi);
% xline(0.1, 'LineWidth', 2);
% xline(0.2, 'LineWidth', 2);
% line([0.1, 0.2], [35, 35], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% x = [0.5 0.68];
% y = [0.8 0.8];
% annotation('textarrow',x,y,'String','Zona de aumento de lesão')
% x = [0.5 0.7];
% y = [0.3 0.6];
% annotation('textarrow',x,y,'String','Sem lesão')
% xlim([x1, x2]);
% ylim([y1, y2]);


% polos = pole(FT);
% zeros = zero(FT);
% format long
% 
% disp('Polos da função de transferência:');
% disp(polos);
% disp('Zeros da função de transferência:');
% disp(zeros);

% Coeficientes do numerador e denominador
Y = [0.02857, 5.526, 7.483e5, 4.914e7, 5.94e9, 1.522e11, 6.092e11, -0.09089, -2.809e-15];
X = [1, 206.8, 2.621e7, 2.044e9, 2.872e11, 1.052e13, 4.817e14, 9.119e15, 3.425e16];

% Definindo a função de transferência
H = tf(Y, X);

% Escolher um zero para o compensador próximo aos polos dominantes
z_c = -5;  % Posição do zero do compensador

% Escolher um polo que mova os polos dominantes para a esquerda
p_c = -9;  % Posição do polo do compensador

% Definir a função de transferência do compensador
Gc = tf([1, -z_c], [1, -p_c]);

% Sistema compensado
H_compensado = series(Gc, H);

% Calcular o ganho necessário para atingir os polos desejados
dominant_pole = -5.73 + 43.01i ;  % exemplo de polo dominante do sistema original
eval_H_compensado = evalfr(H_compensado, dominant_pole);
K = 1 / abs(eval_H_compensado);

% Aplicar o ganho ao sistema compensado
H_compensado = K * H_compensado;

% Verificar os polos de malha fechada resultantes
poles = pole(feedback(H_compensado, 1));

% Exibir o ganho e os polos resultantes
disp('Ganho do compensador:');
disp(K);
disp('Polos de malha fechada resultantes:');
disp(poles);

% Analisar o lugar das raízes do sistema compensado
figure;
rlocus(H_compensado);
title('Lugar das Raízes do Sistema Compensado');

% Resposta ao degrau do sistema original e do sistema compensado
figure;
step(feedback(H, 1));
hold on;
step(feedback(H_compensado, 1));
legend('Sistema Original', 'Sistema Compensado');
title('Resposta ao Degrau');



T = feedback(H, 1);
T2 = feedback(H_compensado, 1);


tempoDeSimulacao = 2; %s
amostra = 0.0001;
tempo = 0:amostra:tempoDeSimulacao;
x0 = [0; 0; 0; 0; 0; 0; 0; 0];

u = 55438.083333333*ones(size(tempo));

[y_pi2, t_pi2] = lsim(T2, u, tempo, x0);
y_pi2 = y_pi2/9.81;

% Plot comparison of time responses
figure(4)
plot(t_pi, y_pi, t_pi2, y_pi2)
ylim([-150 150])
xlim([0 1])
title('Comparação sistema em malha fechada com e sem compensador');
legend('PI - Lugar das Raízes', 'PI com Compensador')
ylabel('Aceleração (Gs)')
xlabel('Tempo(s)')
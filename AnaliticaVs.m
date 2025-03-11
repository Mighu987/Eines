%% Parámetros del problema
clear; clc; close all;

d = 0.08;                  % Distancia entre placas (m)
V0 = 1;                    % Potencial en la placa izquierda (V)
epsilon_r = 1;             % Constante dieléctrica relativa
epsilon_0 = 8.854e-12;     % Permisividad del vacío (F/m)
rho_0 = 1e-8;             % Densidad de carga volumétrica (C/m^3)
num_nodes = 8;            % Número de nodos para la solución FEM

%% Discretización del dominio (FEM)
x = linspace(0, d, num_nodes); % Puntos nodales FEM
h = x(2) - x(1);              % Tamaño de elemento FEM

%% Inicialización de la matriz global y el vector de fuerzas
K = zeros(num_nodes);   % Matriz de rigidez global (FEM)
F = zeros(num_nodes, 1);% Vector de fuerzas global (FEM)

%% Ensamblaje de la matriz FEM
for e = 1:num_nodes-1
    % Matriz elemental FEM
    Ke = (epsilon_r * epsilon_0 / h) * [1, -1; -1, 1];  
    % Vector de fuerzas elemental FEM
    Fe = (-rho_0 * h / 2) * [1; 1];                     
    
    % Ensamblaje en la matriz global FEM
    K(e:e+1, e:e+1) = K(e:e+1, e:e+1) + Ke;
    F(e:e+1) = F(e:e+1) + Fe;
end

%% Aplicación de condiciones de frontera en FEM
K(1, :) = 0; K(1, 1) = 1; F(1) = V0;     % Placa izquierda con potencial fijo
K(end, :) = 0; K(end, end) = 1; F(end) = 0; % Placa derecha con potencial cero

%% Resolución del sistema FEM para obtener el potencial en los nodos
V = K \ F;  

%% Cálculo del campo eléctrico en FEM
E = -diff(V) / h;
x_E = x(1:end-1) + h/2; % Puntos intermedios donde se evalúa el campo FEM  

%% Solución analítica del potencial y el campo eléctrico
x_fine = linspace(0, d, 1000); % Malla fina para la solución analítica
V_exact_fine = (rho_0 / (2 * epsilon_r * epsilon_0)) * x_fine.^2 ...
               - (rho_0 * d / (2 * epsilon_r * epsilon_0) + V0 / d) * x_fine + V0;
E_exact_fine = (rho_0 / (epsilon_r * epsilon_0)) * (d / 2 - x_fine);

%% Gráficos de los resultados
figure;

% Gráfico del potencial eléctrico
subplot(2, 1, 1);
plot(x, V, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); % FEM
hold on;
plot(x_fine, V_exact_fine, 'r--', 'LineWidth', 1.5); % Analítico
legend('FEM', 'Analítica', 'Location', 'best');
xlabel('x (m)'); ylabel('Potencial Eléctrico (V)');
title('Distribución del Potencial Eléctrico');
grid on;

% Gráfico del campo eléctrico
subplot(2, 1, 2);
plot(x_E, E, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); % FEM
hold on;
plot(x_fine, E_exact_fine, 'r--', 'LineWidth', 1.5); % Analítico
legend('FEM', 'Analítica', 'Location', 'best');
xlabel('x (m)'); ylabel('Campo Eléctrico (V/m)');
title('Distribución del Campo Eléctrico');
grid on;

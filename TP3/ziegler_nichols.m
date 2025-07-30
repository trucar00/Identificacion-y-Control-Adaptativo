%Codigo para encontrar K_c y T_c

K_list = 12:0.2:20;
Data = [];
T = [];
sustained_oscillation = false;
Kc = 0; % Ganancia crítica
Tc = 0; % Período de la oscilación

for i = 1:size(K_list, 2)
    K = K_list(i);
    out = sim("practico3_ziegler.slx");
    Data{i} = out.simout.Data; % Salida
    T{i} = out.simout.Time; % Tiempo

    % Encontrar los picos y sus tiempos correspondientes
    [peaks, locs] = findpeaks(Data{i}, T{i});
    rate_list = [];
    
    % Calcular la tasa de cambio y verificar si la oscilación es sostenida
    rate_list = diff(peaks(2:end)) ./ diff(locs(2:end));

    threshold = -5e-05;
    avg_rate = mean(rate_list);
    fprintf('avg= %d\n', avg_rate);
    if avg_rate >= threshold
        sustained_oscillation = true;
        Kc = K;
        Tc = mean(diff(locs(2:end)));
        break
    end
end

% Calcular los parámetros
if sustained_oscillation
    Kp = 0.6*Kc;
    Ki = 1.2*Kc/Tc;
    Kd = 0.075*Kc*Tc;
else
    warning('No lo funcionó.');
end
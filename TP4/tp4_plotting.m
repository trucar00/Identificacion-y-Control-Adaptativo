% Make plots for TP4

y = out.nivout.Data;
t = out.nivout.Time;
ref = ones(1, length(y))*h_list(4);

plot(t, y, t, ref,"--", "LineWidth", 2);
grid on;
xlabel("Tiempo (s)");
ylabel("Nivel del tanque (m)");
legend("Respuesta", "Referencia", "FontSize", 15);
%axis([0 600 0.1 0.9])
title("Controlador no adaptativo");
set(gca, 'FontSize', 15);

% Add a textbox with constants
%dim = [0.6 0.25 0.3 0.3]; % [x y width height] in normalized figure coordinates
%str = sprintf('Parámetros:\nK_p = %.2f\nK_i = %.2f\nK_d = %.2f', Kp, Ki, Kd); % Text to display
%annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
%    'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 15);
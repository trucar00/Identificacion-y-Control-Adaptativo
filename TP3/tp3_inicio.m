%%Initialization
r_0 = 0.05; %Radio menor
r_o = 0.0254; %Radio de orifício de salida
g = 9.81;
theta = 60;

a = pi*r_o^2; %La sección de tubería de salida

F_in = 0.0035;

h_list = [0.2, 0.4, 0.6, 0.8];
k_list = zeros(1, length(h_list));

for i = 1:length(h_list)
    k_list(i) = F_in / (a*sqrt(2*9.81*h_list(i)));
end

F_out_sin_k = a*sqrt(2*g);

tau = 10;
step_fin = 1;

[A, B, C, D] = tf2ss([0 1], [10 1]); % Matrices por las dinámicas
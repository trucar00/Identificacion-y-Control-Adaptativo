clear all;

r_0 = 0.05; %Radio menor
r_o = 0.0254; %Radio de orifício de salida
g = 9.81;
theta = 60;

a = pi*r_o^2; %La sección de tubería de salida

F_in = 0.0035;

%k_out = 1; %Porcentaje de apertura

h_bar = 0.2; %punto de trabajo
k_bar = F_in / (a*sqrt(2*9.81*h_bar)); %k en el punto de trabajo

F_out = a*sqrt(2*g);

%%Linearización
syms h k

F_out_lin = k*a*sqrt(2*g);

f = (F_in - F_out_lin*sqrt(h)) / (pi * ((h^2)/3 + 2*r_0*h/sqrt(3) + r_0^2));

dfdh = diff(f, h);
dfdh_eval = double(vpa(subs(dfdh, [h, k], [h_bar, k_bar])));

dfdk = diff(f, k);
dfdk_eval = double(vpa(subs(dfdk, [h, k], [h_bar, k_bar])));

k_list = [0.05, 0.15, 0.30, 0.50, 0.90];

%%Discretización
Ts = 1; %tiempo del muestreo

A_d = exp(dfdh_eval*Ts);
B_d = dfdh_eval^(-1) * (A_d-1) * dfdk_eval;

%%Tustin
G = tf([0 dfdk_eval], [1 -dfdh_eval]);
Gz = c2d(G, Ts, 'Tustin');
[num,den] = tfdata(Gz, 'v');


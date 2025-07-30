% Separación de dinamicas

h_bar = h_list(1); %punto de trabajo
k_bar = F_in / (a*sqrt(2*g*h_bar)); %k en el punto de trabajo

F_out = a*sqrt(2*g);

%%Linearización
syms h k

F_out_lin = k*a*sqrt(2*g);

f = (F_in - F_out_lin*sqrt(h)) / (pi * ((h^2)/3 + 2*r_0*h/sqrt(3) + r_0^2));

dfdh = diff(f, h);
dfdh_eval = double(vpa(subs(dfdh, [h, k], [h_bar, k_bar])));

dfdk = diff(f, k);
dfdk_eval = double(vpa(subs(dfdk, [h, k], [h_bar, k_bar])));

%%Discretización
Ts = 0.1; %tiempo del muestreo

%% Tustin
G = tf([0 dfdk_eval], [1 -dfdh_eval]);
Gz = c2d(G, Ts, 'Tustin');
[B_d,A_d] = tfdata(Gz, 'v');

%% Ecuacion diofantina

% Resolver por r0 y s0
syms r0 s0   

b0d = B_d(1); b1d = B_d(2);
a1d = A_d(2);

alpha = -0.5; % A_lc(q) = 1+alpha*q^{-1}

eq1 = r0 + b0d * s0 == 1;
eq2 = a1d * r0 + b1d * s0 == alpha;

sol = solve([eq1, eq2], [r0, s0]);

R = double(sol.r0); % R = r0
S = double(sol.s0); % S = s0

disp('R(q) = '), disp(R);
disp('S(q) = '), disp(S);

% T(q)

A_lc = [1 alpha];
T = tf(fliplr(A_lc), fliplr(B_d), Ts);
%disp('T(q) = '), disp(T);

% estimate parameters a1 b0 and b1 with the sinus input, estimate y(k) for each
% iteration with the parameters, r(k)=const, find u(k).
% Step 1: Run RLS with sine-wave input to estimate parameters
% (Do this separately before running STR)

% Step 2: Use these parameters in STR with constant reference
ref = ones(N,1)*h_list(3);
u = zeros(N,1);
y = zeros(N,1);
y(1) = h_list(1); % initial height
u(1) = k_list(1); % initial valve opening

alpha = -0.5; % alc1
t0 = -(1+alpha)/b1;
s0 = -(1/b1)*(alpha-a1);

K_plant = b1/(1-a1);
uss = -0.6/K_plant;
for k = 2:N
    
    y(k) = a1 * y(k-1) + b1 * u(k-1);

    u(k) = t0 * ref(k) - s0 * y(k); % + uss;
    if u(k) < 0
        u(k) = 0;
    elseif u(k) > 1
            u(k) = 1;
    end
    
end

time = (0:N-1) * 10;
u_timeseries = timeseries(u, time);
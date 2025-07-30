    % Method of least squares, direct
    
    % Subtract the initial output and input
    y = y-y(1);
    u = u-u(1);
    
    n = 1; % Output order
    m = 1; % Input order
    
    N = length(y); % Number of samples
    
    % Preallocate Phi and Y
    Phi = zeros(N-1, n + m);
    Y = zeros(N-1, 1);
    
    % Construct Phi and Y matrices
    for k = 2:N
        phi_y = y(k-1);
        phi_u = u(k-1);
        Phi(k-1, :) = [phi_y, phi_u];
        Y(k-1) = y(k);
    end
    
    theta_hat = (Phi' * Phi) \ (Phi' * Y);
    
    parameters = theta_hat;
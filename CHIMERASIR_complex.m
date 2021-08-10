function sol = CHIMERASIR_complex(params, domain, ins)

M = zeros(3, domain(2) + 1);

%Defne initial conditions
M(1,1) = params(1);         
M(2,1) = params(2);
M(3,1) = params(3);
M(4,1) = params(3);
M(5,1) = params(4);

%Define parameters
gamma = params(5);
mu_h = params(6);
z = params(7);
mu_delta = params(8);

%Total population

N = sum(M(:, 1));

for i = domain(1) + 1 : domain(2)    

    %Define probability of interaction
    if M(1, i) > 1        
        red = z * M(2, i) * (M(1, i)/N);
        Probability = ( 1 - 1/M(1, i)) ^ red;
    else        
        Probability = 0;        
    end
    
    % Susceptible
    M(1, i + 1) = M(1, i) * Probability * (1 - mu_h) + N * mu_h;
    % Infected
    M(2, i + 1) = M(1, i) * (1 - Probability) * (1 - mu_h) + M(2, i) * (1 - exp( - gamma)) *  (1-mu_delta) * (1 - mu_h);
    % Recovered
    M(3, i + 1) = M(3, i) * (1 - mu_h) + M(2, i) * exp( - gamma) *(1-mu_delta)* (1 - mu_h);    
    M(4, i + 1) = M(4, i) + M(2, i) * exp( - gamma) * (1 - mu_h)*(1-mu_delta);
    M(5, i + 1) = M(5, i) + M(2, i) *  mu_delta * (1 - mu_h);
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = M; 

end
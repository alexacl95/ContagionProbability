function sol = CHIMERASIR_simple(params, domain, ins)

M = zeros(3, domain(2) + 1);


%Defne initial conditions
M(1,1) = params(1);         
M(2,1) = params(2);
M(3,1) = params(3);
M(4,1) = params(3);
M(5,1) = params(4);

%Define parameters
gamma = params(5);
mu = params(6);
z = params(7);
mu_delta = params(8);

%Total population
N = sum(M(:, 1));

for i = domain(1) + 1 : domain(2)
    
    %Define probability of interaction
    
    
    % Susceptible
    M(1, i + 1) = M(1, i) * (1 - mu) * psi + mu * N;
    % Infected
    M(2, i + 1) = M(1, i) * (1 - mu) * (1 - psi) + M(2, i) * (1 - exp( - gamma)) *  (1-mu_delta) * (1 - mu);
    % Recovered
    M(3, i + 1) = M(3, i) * (1 - mu) + M(2, i) * exp( - gamma) * (1 - mu);
    M(4, i + 1) = M(4, i) + M(2, i) * exp( - gamma) * (1 - mu)*(1-mu_delta);
    M(5, i + 1) = M(5, i) + M(2, i) *  mu_delta * (1 - mu);
    
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = M; 

end
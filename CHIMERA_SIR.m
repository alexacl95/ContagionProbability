function sol = CHIMERA_SIR(params, domain, ins)

M = zeros(4, domain(2) + 1);

%Defne initial conditions
M(1,1) = params(1);         
M(2,1) = params(2);
M(3,1) = params(3);
M(4,1) = params(3);
M(5,1) = params(4);
M(6,1) = sum(M(1:3,1));

%Define parameters
gamma = params(5);
mu_h = params(6);
z = params(7);
nu = params(8);
%Total population

prob = ins.Prob;

for i = domain(1) + 1 : domain(2)    
    
    N = sum(M(1:3,i));
    
    aleph = 1 + nu*M(2,i)/(M(1,i)+ M(2,i));
    
    %Define probability of interaction
    switch prob 
        case 1 %psi
            Probability = (1 - ((M(2, i)/N))^aleph)^z;
        case 2 %phi
                red = (z * M(2, i) * M(1, i)/N);
                Probability = (1-1/(M(1, i)+1)) ^ red^aleph;
        case 3 %classic
            Probability = 1 - z*(M(2,i)/N)^aleph;
    end
    
    % Susceptible
    M(1, i + 1) = M(1, i) * Probability + M(3,i) * mu_h;
    
    % Infected
    M(2, i + 1) = M(1, i) * (1 - Probability) +...
        M(2, i) * (1 - gamma);
    
    % Recovered
    M(3, i + 1) = M(3, i) * (1 - mu_h) + M(2, i) * gamma;    
    
    %Fitting functions
    %%accumulated recovered
    M(4, i + 1) = M(4, i) + M(2, i) *gamma ;
    %%accumulated cases
    M(5, i + 1) = M(5, i) + M(1, i) * (1 - Probability);
    %%Total popultion
    M(6, i + 1) = sum(M(1:3,i));
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = M; 

end
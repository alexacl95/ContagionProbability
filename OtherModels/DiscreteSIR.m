function sol = DiscreteSIR(params, domain, ins)

M = zeros(4, domain(2) + 1);

%Defne initial conditions
M(1,1) = params(1);         
M(2,1) = params(2);
M(3,1) = params(3);
M(4,1) = params(3);
M(5,1) = sum(M(1:3,1));

%Define parameters
gamma = params(4);
alpha = params(5);
A = params(6); 
B = params(7);
delta = params(8);
sigma = params(9);


prob = ins.Prob;

for i = domain(1) + 1 : domain(2)    

    %Define probability of interaction
    
    N = sum(M(1:3,i));
    x = alpha*M(2,i)/N;
    
    switch prob 
        case 1
            Probability = exp(-x);
        case 2
            Probability = B/(B+x);
    end
    
    A = (1-gamma)*(sum(M(1:3,i)));
    
    % Susceptible
    M(1, i + 1) = A + gamma*M(1, i)*Probability + gamma*(1-delta)*M(3, i);
    
    % Infected
    M(2, i + 1) = gamma*M(1, i)*(1-Probability) + gamma*sigma*M(2,i);
    
    % Recovered
    M(3, i + 1) = gamma*(1-sigma)*M(2,i) + gamma*delta*M(3,i);
    
    %Fitting functions
    %%accumulated recovered
    M(4, i + 1) = M(4, i) + gamma*(1-sigma)*M(2,i);
    M(5, i+1) = sum(M(1:3,i));
end

    sol = struct();
    sol.x = domain(1) : domain(2);
    sol.y = M; 

end
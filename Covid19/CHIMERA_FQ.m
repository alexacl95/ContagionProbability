function sol = CHIMERA_FQ(params, domain, ins)

M = zeros(11, domain(2) + 1);
% Defne initial conditions

M(1,1) = params(1);         % S_f        
M(2,1) = 0;                 % S_q
M(3,1) = params(2);         % I_f
M(4,1) = params(3);         % I_q
M(5,1) = params(4);         % I_j
M(6,1) = params(5);         % R_total
M(7,1) = params(6);         % R_j
M(8, 1) = sum(M(1:7,1));    % N
M(11,1) = params(4);        % Accumulated cases    


% Define parameters
lambda = params(7); % Enter quarantine
alpha = params(8);  % Leave quarantine
mu = params(9);    % Leave R
theta = params(10); % Detection
gamma = params(11); % Recovery
nu = params(12);    % conectivity    
z = params(13);     % Interactions
beta =  params(14);

prob = ins.Prob;

for i = domain(1) + 1 : domain(2)    
    
    N =  sum(M(1:6, i));
    I_tot = sum(M(3:5,i));
        
    %% Probability of interaction
    
    aleph = 1 + nu*I_tot/(M(1,i) + I_tot);
    
    switch prob 
        case 1 %psi
            Probability = beta*(1-(1-(M(3,i)/N)^aleph)^z);
        case 2 %phi
            if M(1, i)>=1
                Probability = beta*(1-(1-1/(M(1, i)))^(z *M(3,i)*(M(1, i)/N)^aleph));
            else
                Probability=0;
            end

        case 3 %classic
            Probability = beta*(M(3,i)/N)^aleph;
        case 4
            Probability = beta*(M(3,i)/N);
    end
    %% 
    % Susceptible
      M(1, i + 1) = (1 - lambda)*(1-Probability)*M(1, i) + alpha*M(2, i) + mu*M(6, i);
      M(2, i + 1) = lambda*(1-Probability)*M(1,i) + (1-alpha)*M(2, i);
      
    % Infected
      M(3, i + 1) = Probability*M(1,i) + (1-lambda)*(1-theta)*(1-gamma)*M(3,i) +...
          alpha*(1 - theta)*(1-gamma)*M(4,i);
      M(4, i + 1) = lambda*(1 - gamma)*(1 - theta)*M(3, i) + (1 - gamma)*(1-theta)*(1-alpha)*M(4,i);
      M(5, i + 1) = theta*(1-gamma)*M(3, i) + theta*(1 - gamma)*M(4, i) + M(5,i)*(1-gamma);
      
    % Recovered
      M(6, i + 1) = M(6, i)*(1 - mu) + gamma*(M(3, i) + M(4, i) + M(5, i));
    
      
    %% Secondary functions   
    % Recovered identified    
      M(7, i + 1) = M(7, i)*(1 - mu) + M(5,i)*gamma;
      
    % Accumulated cases
      M(11, i + 1) = M(11, i) + theta*(1-gamma)*M(3, i) + theta*(1 - gamma)*M(4, i);
      
    % Total popultion
      M(8, i + 1) = N;
      
    % Probabilities  
      M(9, i + 1) = aleph;
      M(10,i + 1) = Probability;
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = M; 

end
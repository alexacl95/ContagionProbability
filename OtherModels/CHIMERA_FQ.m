function sol = CHIMERA_FQ(params, domain, ins)

M = zeros(11, domain(2) + 1);

%Defne initial conditions
M(1,1) = params(1);         
M(2,1) = params(2);
M(3,1) = params(3);
M(4,1) = params(4);
M(5,1) = params(5);
M(6,1) = params(6);
M(7,1) = params(7);
M(11,1) = params(4);

%Define parameters
lambda = params(8); %Enter quarantine
alpha = params(9);  %Leave quarantine
mu = params(10);    %Leave R
theta = params(11); %Detection
gamma = params(12); %Recovery
nu = params(13);
z = params(14);

prob = ins.Prob;

M(8, 1) = sum(M(1:7,1));

for i = domain(1) + 1 : domain(2)    
    
    N =  sum(M(1:7, i));
%     I_tot = sum(3:5,i);
    
    aleph = 1 + nu*M(3,i)/(M(1,i)+ M(3,i));
%     aleph = 1 + nu*I_tot/(M(1,i) + I_tot);
%       aleph = 1 + nu*M(3,i)/(M(1,i)/nu+ M(3,i));
    
    %Define probability of interaction
    switch prob 
        case 1 %psi
            Probability = (1 - ((M(3, i)/N))^aleph)^z;
        case 2 %phi
            if M(1, i) > 1        
                red = (z * M(3, i) * M(1, i)/N);
                Probability = ((M(1, i) - 1)/M(1, i)) ^ red^aleph;
            else        
                Probability = 0;        
            end
        case 3 %classic
            Probability = 1 - z*(M(3,i)/N)^aleph;
    end
    
    % Susceptible
      M(1, i + 1) = (1 - lambda)*Probability*M(1, i) + alpha*M(2, i) + mu*(M(6, i) + M(7, i));
      M(2, i + 1) = lambda*Probability*M(1,i) + (1-alpha)*M(2, i);
      
    % Infected
      M(3, i + 1) = (1 - Probability)*M(1,i) + (1-lambda)*(1-theta)*(1-gamma)*M(3,i) +...
          alpha*(1 - theta)*(1-gamma)*M(4,i);
      M(4, i + 1) = lambda*(1 - gamma)*(1 - theta)*M(3, i) + (1 - gamma)*(1-theta)*(1-alpha)*M(4,i);
      M(5, i + 1) = theta*(1-gamma)*M(3, i) + theta*(1 - gamma)*M(4, i) + M(5,i)*(1-gamma);
      
    % Recovered
      M(6, i + 1) = M(6, i)*(1 - mu) + gamma*(M(3, i) + M(4, i));
      M(7, i + 1) = M(7, i)*(1 - mu) + M(5,i)*gamma;
    
    %Fitting functions
    %%accumulated recovered
   
    %%accumulated cases
    M(11, i + 1) = M(11, i) + theta*(1-gamma)*M(3, i) + theta*(1 - gamma)*M(4, i);
    %%Total popultion
    M(8, i + 1) = sum(M(1:7, i));
    M(9, i + 1) = aleph;
    M(10,i + 1) = Probability;
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = M; 
% % M(8,end)
end
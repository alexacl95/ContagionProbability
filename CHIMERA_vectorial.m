function sol = CHIMERA_vectorial(params, domain, ins)

domain(1) = 0;

M = zeros(6, domain(2) + 1);

%Defne initial conditions
M(1, 1) = params(1);         %Hs
M(2, 1) = params(2);         %Hi
M(3, 1) = params(3);         %Hr
M(4, 1) = params(4);         %Ms
M(5, 1) = params(5);         %Mi
M(6, 1) = params(2);         %Hi_inst

%Define parameters

gamma = params(7);
mu_m = params(8);
mu_h = params(9);
z_m = params(10);
z_h = params(11);
r = params(12);
C = params(13);

%Total population
N_M = zeros(1, domain(2) + 1);

for i = domain(1) + 1 : domain(2)
    
    %Total Population H
    N_H = sum(M(1 : 3, i));
    
    %Define probability of interaction
    if M(1, i) > 1        
        %to find a mosquito;
        red_m = z_m * M(5, i) * (M(1, i)/N_H);
        P_m = ((M(1, i) - 1)/M(1, i)) ^ red_m;             
    else        
        P_m = 0;  
    end
    
    %to find a human
    P_h = ((N_H - M(2, i))/N_H) ^ z_h;
    
    %% Humans
    % Susceptible
    M(1, i + 1) = M(1, i) * P_m * (1 - mu_h) + (M(1, i) +  M(2, i) + M(3, i)) * mu_h;
    % Infected
    M(2, i + 1) = M(1, i) * (1 - P_m) * (1 - mu_h) + M(2, i) * (1 - exp( - gamma)) * (1 - mu_h);
    M(6, i + 1) = M(1, i) * (1 - P_m) * (1 - mu_h);
    % Recovered
    M(3, i + 1) = M(3, i) * (1 - mu_h) + M(2, i) * exp( - gamma) * (1 - mu_h);
    
    %% Mosquitoes
    % Recruitment rate
    
    %Total population
    N_M(i) = sum(M(4 : 5, i));
    
    if i > 2
        A = r * N_M(i - 2) * exp( (1 - N_M(i)/C)); %Ricker
    else
        A = r * N_M(1) * exp( (1 - N_M(1)/C)); %Ricker
    end
    
    % Susceptible
    M(4, i + 1) = M(4, i) * P_h * (1 - exp( - mu_m)) + A;
    % Infected
    M(5, i + 1) = M(4, i) * (1 - P_h) * (1 - exp( - mu_m)) + M(5, i) * (1 - exp(- mu_m));
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = [M; cumsum(M(6,:))]; 

end
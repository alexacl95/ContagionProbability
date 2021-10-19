function sol = CHIMERA_Vector(params, domain, ~)

M = zeros(8, domain(2) + 1);

%Defne initial conditions
M(1, 1) = params(1);         %Hs
M(2, 1) = params(2);         %Hi
M(3, 1) = params(3);         %Hr
M(4, 1) = params(4);         %Ms
M(5, 1) = params(5);         %Mi
M(6, 1) = params(6);         %Md %quitar esta
M(7, 1) = params(2);         %Hi_inst
M(8, 1) = sum(M(1:3,1));
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
        phi = ((M(1, i) - 1)/M(1, i)) ^ red_m;             
    else        
        phi = 0;  
    end
    
    %to find a human
    psi = ((N_H - M(2, i))/N_H) ^ z_h;
    
    %% Humans
    % Susceptible
    M(1, i + 1) = M(1, i) * phi * (1 - mu_h) + (M(1, i) +  M(2, i) + M(3, i)) * mu_h;
    % Infected
    M(2, i + 1) = M(1, i) * (1 - phi) * (1 - mu_h) + M(2, i) * (1 - gamma) * (1 - mu_h);
    % Recovered
    M(3, i + 1) = M(3, i) * (1 - mu_h) + M(2, i) * gamma * (1 - mu_h);
    
    %fitting curve: instantaneous infections
    M(7, i + 1) = M(1, i) * (1 - phi) * (1 - mu_h);
    
    %% Mosquitoes
    %Total population
    N_M(i) = sum(M(4 : 5, i));    
    
    % Recruitment rate
    if i > 2
        A = r * N_M(i - 2) * exp( (1 - N_M(i)/C)); %Ricker
    else
        A = r * N_M(1) * exp( (1 - N_M(1)/C));
    end
    
    % Susceptible
    M(4, i + 1) = M(4, i) * psi * (1 - mu_m) + A;
    % Infected
    M(5, i + 1) = M(4, i) * (1 - psi) * (1 - mu_m) +...
        M(5, i) * (1 - mu_m);
    % Dead
    M(6, i + 1) = M(5, i) * mu_m +...
        M(4, i) * mu_m; %Not accumulated
    
    M(8, i + 1) = N_H;
end

sol = struct();
sol.x = domain(1) : domain(2);
sol.y = [M; cumsum(M(7,:))]; 

end
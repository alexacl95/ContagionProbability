function sol = CHIMERA_FQ_continuous(params, domain, ins)

IC = [    
     params(1);  %I_j    
     params(2);  %R_j
     params(3);  %S_f
     params(4);  %S_q    
     params(5);  %I_f                
     params(6);  %I_q
     params(7);  %R
    ];

p = [params(8);     %alpha
     params(9);     %beta
     params(10);    %gamma 
     params(11);    %lambda 
     params(12);    %mu 
     params(13);    %nu             
     params(14);    %theta
     params(15);    %z  
     ins.Prob
];

% options = odeset('NonNegative',1:7,'RelTol',1e-5,'AbsTol',1e-5,'MaxStep',1e-2);

[t,y] = ode45(@(t,y) f(t,y, p), [domain(1) domain(2)], IC);

N = sum(y(1,:));
z = params(15);
I_f = y(:,3);
S_f = y(:,1);
aleph = 1 + params(13)*I_f./(I_f + S_f);


    switch ins.Prob 
        case 1 %psi
            P = 1-(1 - (I_f/N).^aleph).^z;
        case 2 %phi
            P = 1-(1 - 1./(S_f+1)).^((z*I_f).*(S_f./N).^aleph);
        case 3 %classic
            P = (I_f/N).^aleph;
    end 

sol = struct();
sol.x = t;
sol.y = [y, aleph, P]'; 

end

function dX = f(~, X, params)
    S_f = X(1);
    S_q = X(2);
    I_f = X(3);
    I_q = X(4);
    I_j = X(5);
    R = X(6);
    R_j = X(7);
    
    lambda = params(1); %Enter quarantine
    alpha = params(2);  %Leave quarantine
    mu = params(3);    %Leave R
    theta = params(4); %Detection
    gamma = params(5); %Recovery
    nu = params(6);
    z = params(7);
    beta = params(8);
    prob = params(9);
    
    N = S_f + S_q + I_f + I_q + I_j + R + R_j;
    aleph = 1 + nu*I_f/(I_f + S_f);
    
    switch prob 
        case 1 %psi
            P = 1-(1 - (I_f/N)^aleph)^z;
        case 2 %phi
            P = 1-(1 - 1/(S_f+1))^((z*I_f)*(S_f/N)^aleph);
        case 3 %classic
            P = (I_f/N)^aleph;
    end

    dX = [-S_f*(beta*P + lambda) + S_q*alpha + mu*(R_j + R);
          lambda*S_f - S_q*alpha;
          S_f*beta*P + alpha* I_q - I_f*(gamma + theta + lambda);
          -I_q*(alpha + gamma + theta) + I_f*lambda;
          theta*(I_q + I_f) - gamma*I_j;
          (I_q + I_f)*gamma - mu*R;
          I_j*gamma - mu*R_j
          ];

end


function dx = plumage_preference(t,x,paramVec)    
    
    % memory allocation
    dx = zeros(8,1);
    
    M_AB = x(1);
    M_Ab = x(2);
    M_aB = x(3);
    M_ab = x(4);
    F_AB = x(5);
    F_Ab = x(6);
    F_aB = x(7);
    F_ab = x(8);
    
    M = sum(x(1:4));
    F = sum(x(5:end)); 
    N = M + F; 
    
    % model parameters
    % alpha parameters describe extent of preference 
    alpha_A = paramVec(1); 
    alpha_a = paramVec(2);
    lambda = paramVec(3);
    K = paramVec(4);
    mu = paramVec(5);
    fitness_M_AB = paramVec(6);
    fitness_M_Ab = paramVec(7);
    fitness_M_aB = paramVec(8);
    fitness_M_ab = paramVec(9);
    fitness_F_AB = paramVec(10);
    fitness_F_Ab = paramVec(11);
    fitness_F_aB = paramVec(12);
    fitness_F_ab = paramVec(13);
    
    chi_A = alpha_A*(M_ab + M_aB)/(M_AB+M_Ab); % unfavorable/favorable pairings
    chi_a = alpha_a*(M_Ab + M_AB)/(M_aB+M_ab);
    
    % genotype probabilities
    p_AB = (1/(M*F))*((1+chi_A)*(M_AB*F_AB + (1/2)*M_AB*F_aB + (1/2)*M_Ab*F_AB + (1/4)*M_Ab*F_aB) + ...
        (1+chi_a)*((1/4)*M_aB*F_Ab) + ...
        (1-alpha_A)*((1/2)*M_aB*F_AB + (1/4)*M_ab*F_AB) + ...
        (1-alpha_a)*((1/2)*M_AB*F_Ab + (1/4)*M_AB*F_ab));
    p_Ab = (1/(M*F))*((1+chi_A)*((1/2)*M_Ab*F_AB + (1/4)*M_Ab*F_aB) + ...
        (1+chi_a)*((1/4)*M_aB*F_Ab + (1/2)*M_ab*F_Ab) + ...
        (1-alpha_a)*((1/4)*M_AB*F_ab + M_Ab*F_Ab + (1/2)*M_AB*F_Ab + (1/2)*M_Ab*F_ab) + ...
        (1-alpha_A)*((1/4)*M_ab*F_AB));
    p_aB = (1/(M*F))*((1+chi_A)*((1/2)*M_AB*F_aB + (1/4)*M_Ab*F_aB) + ...
        (1+chi_a)*((1/4)*M_aB*F_Ab + (1/2)*M_aB*F_ab) + ...
        (1-alpha_a)*((1/4)*M_AB*F_ab) + ...
        (1-alpha_A)*((1/2)*M_aB*F_AB + M_aB*F_aB + (1/4)*M_ab*F_AB + (1/2)*M_ab*F_aB));
    p_ab = (1/(M*F))*((1+chi_A)*((1/4)*M_Ab*F_aB) + ...
        (1+chi_a)*((1/4)*M_aB*F_Ab + (1/2)*M_aB*F_ab + (1/2)*M_ab*F_Ab + M_ab*F_ab) + ...
        (1-alpha_a)*((1/4)*M_AB*F_ab + (1/2)*M_Ab*F_ab) + ...
        (1-alpha_A)*((1/4)*M_ab*F_AB + (1/2)*M_ab*F_aB));
    
    dx(1,1) = lambda*(1/2)*F*p_AB*max((1-N/K),0) - mu*(1+fitness_M_AB)*M_AB; 
    dx(2,1) = lambda*(1/2)*F*p_Ab*max((1-N/K),0) - mu*(1+fitness_M_Ab)*M_Ab; 
    dx(3,1) = lambda*(1/2)*F*p_aB*max((1-N/K),0) - mu*(1+fitness_M_aB)*M_aB; 
    dx(4,1) = lambda*(1/2)*F*p_ab*max((1-N/K),0) - mu*(1+fitness_M_ab)*M_ab;  
    dx(5,1) = lambda*(1/2)*F*p_AB*max((1-N/K),0) - mu*(1+fitness_F_AB)*F_AB;
    dx(6,1) = lambda*(1/2)*F*p_Ab*max((1-N/K),0) - mu*(1+fitness_F_Ab)*F_Ab;
    dx(7,1) = lambda*(1/2)*F*p_aB*max((1-N/K),0) - mu*(1+fitness_F_aB)*F_aB;
    dx(8,1) = lambda*(1/2)*F*p_ab*max((1-N/K),0) - mu*(1+fitness_F_ab)*F_ab;
end
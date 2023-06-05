close all

%% Sandbox section to test different parameter combinations, etc.
alpha_A = 1.0;              % B female preference for A plumage
alpha_a = 0;                % b female preference for a plumage
lambda = 1/3;               % no. of eggs per month (assuming 4 eggs per 
                            %   clutch, 1 clutch per year)
K = 1;                      % carrying capacity
mu = 1/36;                  % (monthly) mortality rate 
maxTime = 30000;            

N = 100;                    % total pop

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
fitnessVec = [0, 0, 0, 0, 0, 0, 0, 0];
init_conditions = [0,1/N,0,(N-2)/(2*N),0,0,1/N,(N-2)/(2*N)];
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
plot(A_male_freq, '-r', 'linewidth', 1.5);
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
set(gca,'fontsize',16);


%% 
% A simple illustrative example of how ornamental plumage evolves within a
% population.
alpha_A = 1.0;              % B female preference for A plumage
alpha_a = 0.0;              % b female preference for a plumage
lambda = 1/3;               % no. of eggs per month (assuming 4 eggs per 
                            %   clutch, 1 clutch per year)
K = 1;                      % carrying capacity
mu = 1/36;                  % (monthly) mortality rate 
maxTime = 30000;            

N = 100;                    % total pop is 2N

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
fitnessVec = [0.1, 0.1, 0, 0, 0, 0, 0, 0];
init_conditions = [0,1/(N),0,(N/2-1)/(N),0,0,1/(N),(N/2-1)/(N)];
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
relative_avg_lifespan = sum((1./(mu*(1+fitnessVec(1:4)))).*x(:,1:4),2) ./ (sum((1/mu)*x(:,1:4),2));

subplot(1,2,1)
plot(A_male_freq, '-r', 'linewidth', 1.5);
title('a.');
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

subplot(1,2,2)
plot(relative_avg_lifespan, '-k', 'linewidth', 1.5);
title('b.');
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('relative avg. male lifespan','interpreter','latex');
set(gca,'fontsize',16);

%% 
% A simple illustrative example of how ornamental plumage can evolve within
% a population, now assuming that there is a slight preference for cryptic
% males when A and B alleles appear, and also assuming preference strength
% for B plumage is reduced 
alpha_A = 1.0;              % B female preference for A plumage
alpha_a = 0.1;                % b female preference for a plumage
lambda = 1/3;               % no. of eggs per month (assuming 4 eggs per 
                            %   clutch, 1 clutch per year)
K = 1;                      % carrying capacity
mu = 1/36;                  % (monthly) mortality rate 
maxTime = 30000;            

N = 100;                    % total pop is 2N

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
fitnessVec = [0.1, 0.1, 0, 0, 0, 0, 0, 0];
init_conditions = [0,1/(N),0,(N/2-1)/(N),0,0,1/(N),(N/2-1)/(N)];
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);

subplot(1,2,1)
plot(A_male_freq, '-r', 'linewidth', 1.5);
title('a.')
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

% second experiment
alpha_A = 0.9;              % B female preference for A plumage
alpha_a = 0.0;              % b female preference for a plumage

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);

subplot(1,2,2)
plot(A_male_freq, '-r', 'linewidth', 1.5);
title('b.');
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

%% 
% How can an ornamental allele become fixed in the male population? 
% As above, but now fitnesses are changed so that the invading allele is
% actually beneficial. 
alpha_A = 1.0;              % B female preference for A plumage
alpha_a = 0.1;                % b female preference for a plumage
lambda = 1/3;               % no. of eggs per month (assuming 4 eggs per 
                            %   clutch, 1 clutch per year)
K = 1;                      % carrying capacity
mu = 1/36;                  % (monthly) mortality rate 
maxTime = 30000;            

N = 100;                    % total pop is 2N

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
fitnessVec = [-0.1, -0.1, 0, 0, 0, 0, 0, 0];
init_conditions = [0,1/(N),0,(N/2-1)/(N),0,0,1/(N),(N/2-1)/(N)];
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
relative_avg_lifespan = sum((1./(mu*(1+fitnessVec(1:4)))).*x(:,1:4),2) ./ (sum((1/mu)*x(:,1:4),2));

subplot(1,2,1)
plot(A_male_freq, '-r', 'linewidth', 1.5);
title('a.')
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

% second experiment
alpha_A = 0.0;              % B female preference for A plumage
alpha_a = 0.0;              % b female preference for a plumage

% fitness costs for each sex and genotype
% M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);

subplot(1,2,2)
plot(A_male_freq, '-r', 'linewidth', 1.5);
title('b.');
hold on
plot(B_female_freq, '-b', 'linewidth', 1.5);
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

%% (Experiments are numbered as they are listed on Page 3 of Notebook 3.)
% Mating preference and fitness cost are varied; plumage and preference
% alleles are recorded for each combination of fitness cost and preference
% strength. 

% Initial conditions are those of a novel mutation, i.e. a single male 
% bird with bright plumage (A) and mating preference for bright plumage (B).
A_fitness_vec = 0.0:0.01:0.5;
alpha_A_vec = 0:0.01:1;

A_freq_male_mat = zeros(length(A_fitness_vec), length(alpha_A_vec)); 
B_freq_female_mat = zeros(length(A_fitness_vec), length(alpha_A_vec)); 

for i = 1:length(A_fitness_vec)
    for j = 1:length(alpha_A_vec)
        disp([i,j]);
        fitnessVec = [A_fitness_vec(i), A_fitness_vec(i), ...
            0, 0, 0, 0, 0, 0];
        % M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
        paramVec = [alpha_A_vec(j), 0, lambda, K, mu, fitnessVec];
        init_conditions = [0,1/N,0,(N-2)/(2*N),0,0,1/N,(N-2)/(2*N)];
        maxTime = 5000;

        [~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

        A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
        B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
        equilBool = range(A_male_freq((end-1000):end));
        equilBool = (equilBool > 0.001);

        counter = 3;
        while (equilBool)
            maxTime = (10^(counter + 1))*5;

            [~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

            A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
            B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
            equilBool = range(A_male_freq((end-1000):end));
            equilBool = (equilBool > 0.001);

            counter = counter + 1; 
        end

        % store results
        A_freq_male_mat(i,j) = A_male_freq(end);
        B_freq_female_mat(i,j) = B_female_freq(end); 
    end
end

% and the plots
[X,Y] = meshgrid(alpha_A_vec,A_fitness_vec);
surf(X,Y,A_freq_male_mat,'FaceColor','interp',...
    'EdgeColor','none');
view(0,90)
ylabel('fitness cost of ornament','interpreter','latex');
xlabel('strength of mating preference, $\alpha_A$','interpreter','latex');
set(gca,'fontsize',16);
colorbar

figure
[X,Y] = meshgrid(alpha_A_vec,A_fitness_vec);
surf(X,Y,B_freq_female_mat,'FaceColor','interp',...
    'EdgeColor','none');
view(0,90)
ylabel('fitness cost of ornament','interpreter','latex');
xlabel('strength of mating preference, $\alpha_A$','interpreter','latex');
colorbar
set(gca,'fontsize',16);

%% Mating preference and fitness cost are varied; now plumage affects
% fitness of females. 

% Initial conditions are those of a novel mutation, i.e. a single male 
% bird with bright plumage (A) and mating preference for bright plumage (B).
A_freq_male_mat_2 = zeros(length(A_fitness_vec), length(alpha_A_vec)); 
B_freq_female_mat_2 = zeros(length(A_fitness_vec), length(alpha_A_vec)); 

for i = 1:length(A_fitness_vec)
    for j = 1:length(alpha_A_vec)
        disp([i,j]);
        fitnessVec = [A_fitness_vec(i), A_fitness_vec(i), ...
            0, 0, A_fitness_vec(i)/2, A_fitness_vec(i)/2, 0, 0];
        % M_AB, M_Ab, M_aB, M_ab, F_AB, F_Ab, F_aB, F_ab
        paramVec = [alpha_A_vec(j), 0, lambda, K, mu, fitnessVec];
        init_conditions = [0,1/N,0,(N-2)/(2*N),0,0,1/N,(N-2)/(2*N)];
        maxTime = 5000;

        [~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

        A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
        B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
        equilBool = range(A_male_freq((end-1000):end));
        equilBool = (equilBool > 0.001);

        counter = 3;
        while (equilBool)
            maxTime = (10^(counter + 1))*5;

            [~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

            A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
            B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
            equilBool = range(A_male_freq((end-1000):end));
            equilBool = (equilBool > 0.001);

            counter = counter + 1; 
        end

        % store results
        A_freq_male_mat_2(i,j) = A_male_freq(end);
        B_freq_female_mat_2(i,j) = B_female_freq(end); 
    end
end

% and the plot
[X,Y] = meshgrid(alpha_A_vec,A_fitness_vec);
surf(X,Y,A_freq_male_mat_2,'FaceColor','interp',...
    'EdgeColor','none');
view(0,90)
ylabel('fitness cost of ornament','interpreter','latex');
xlabel('strength of mating preference, $\alpha_A$','interpreter','latex');
set(gca,'fontsize',16);
colorbar

figure
[X,Y] = meshgrid(alpha_A_vec,A_fitness_vec);
surf(X,Y,B_freq_female_mat_2,'FaceColor','interp',...
    'EdgeColor','none');
view(0,90)
ylabel('fitness cost of ornament','interpreter','latex');
xlabel('strength of mating preference, $\alpha_A$','interpreter','latex');
colorbar
set(gca,'fontsize',16);

%% Does female fitness cost make a difference?
fitnessConst = 0.05;
fitnessVec = [fitnessConst, fitnessConst, 0, 0, 0, 0, 0, 0];
init_conditions = [0,1/N,0,(N-2)/(2*N),0,0,1/N,(N-2)/(2*N)];
maxTime = 20000;
paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
[~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);

subplot(1,2,1)
plot(A_male_freq, 'linewidth', 1.5);
hold on
ylim([0,1]); 
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

subplot(1,2,2)
plot(B_female_freq, 'linewidth', 1.5);
hold on
xlabel('time (months)','interpreter','latex');
ylabel('allele frequencies','interpreter','latex');
set(gca,'fontsize',16);

fracVec = 0.2:0.2:1; 

for i = fracVec
    fitnessVec = [fitnessConst, fitnessConst, 0, 0, fitnessConst*i, fitnessConst*i, 0, 0];
    paramVec = [alpha_A, alpha_a, lambda, K, mu, fitnessVec];
    [~,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);

    A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
    B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);

    subplot(1,2,1)
    plot(A_male_freq, 'linewidth', 1.5);

    subplot(1,2,2)
    plot(B_female_freq, 'linewidth', 1.5);
end

subplot(1,2,1)

subplot(1,2,2)
legend('$f_{F,AB} = f_{F,Ab} = 0.0$','$f_{F,AB} = f_{F,Ab} = 0.01$','$f_{F,AB} = f_{F,Ab} = 0.02$',...
    '$f_{F,AB} = f_{F,Ab} = 0.03$','$f_{F,AB} = f_{F,Ab} = 0.04$','$f_{F,AB} = f_{F,Ab} = 0.05$','interpreter','latex');


%% code for simple plotting of results
simplePlots = false;
if (simplePlots)
    opts = odeset('RelTol',3e-14,'AbsTol',3e-14);
    [t,x]=ode45(@plumage_preference,0:maxTime,init_conditions,opts, paramVec);
    % genotype plots
    plot(t,x(:,1),'--b','LineWidth',1.5);
    hold on
    plot(t,x(:,2),'--g','LineWidth',1.5);
    plot(t,x(:,3),'--k','LineWidth',1.5);
    plot(t,x(:,4),'--m','LineWidth',1.5);
    plot(t,x(:,5),'-b','LineWidth',1.5);
    plot(t,x(:,6),'-g','LineWidth',1.5);
    plot(t,x(:,7),'-k','LineWidth',1.5);
    plot(t,x(:,8),'-m','LineWidth',1.5);
    legend({'$M_{AB}$','$M_{Ab}$','$M_{aB}$','$M_{ab}$','$F_{AB}$',...
        '$F_{Ab}$','$F_{aB}$','$F_{ab}$'},'interpreter','latex','location','best');

    % frequency of A allele in males and total population
    figure
    subplot(2,2,1)
    A_male_freq = (x(:,1) + x(:,2))./ sum(x(:,1:4),2);
    plot(t,A_male_freq,'-r','LineWidth',1.5);
    xlim([0,max(t)]);
    ylim([0,1]);
    title('A freq. in males');

    subplot(2,2,2)
    A_total_freq = (x(:,1) + x(:,2) + x(:,5) + x(:,6))./ sum(x(:,1:8),2);
    plot(t,A_total_freq,'--r','LineWidth',1.5);
    xlim([0,max(t)]);
    ylim([0,1]);
    title('A freq. in total');

    subplot(2,2,3)
    B_female_freq = (x(:,5) + x(:,7))./ sum(x(:,5:8),2);
    plot(t,B_female_freq,'-b','LineWidth',1.5);
    xlim([0,max(t)]);
    ylim([0,1]);
    title('B freq. in females');

    subplot(2,2,4)
    B_total_freq = (x(:,5) + x(:,7) + x(:,1) + x(:,3))./ sum(x(:,1:8),2);
    plot(t,B_total_freq,'--r','LineWidth',1.5);
    ylim([0,1]);
    xlim([0,max(t)]);
    title('B freq. in total');
end


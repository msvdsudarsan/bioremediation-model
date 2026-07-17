clear; clc; close all;

sigma = 1.0;
delta = 0.1;
alpha = 0.8;
K     = 1.0;
theta = 0.6;
mu    = 0.2;
gamma = 0.05;

fprintf('============================================================\n');
fprintf(' MODEL PARAMETERS\n');
fprintf('------------------------------------------------------------\n');
fprintf(' sigma=%.3f  delta=%.3f  alpha=%.3f  K=%.3f\n', sigma,delta,alpha,K);
fprintf(' theta=%.3f  mu=%.3f  gamma=%.3f\n\n', theta,mu,gamma);

J = @(P,M)[ -delta - alpha*M*K/(K+P)^2 , -alpha*P/(K+P);
             M*(theta*K/(K+P)^2 - gamma), theta*P/(K+P) - mu - gamma*P ];

fprintf('============================================================\n');
fprintf(' EQUILIBRIA AND STABILITY\n');
fprintf('------------------------------------------------------------\n');

P0 = sigma/delta;
ev0 = eig(J(P0,0));
fprintf(' Washout  E0: P0=%.4f, M0=0.0000\n', P0);
fprintf('    eigenvalues = [% .4f , % .4f]\n', real(ev0(1)), real(ev0(2)));
fprintf('    --> %s\n\n', classify(ev0));

coef   = [gamma, mu+gamma*K-theta, mu*K];
Proots = roots(coef);
Pstar_list = []; Mstar_list = [];
for i = 1:numel(Proots)
    P = Proots(i);
    if isreal(P) && P>0 && P < sigma/delta
        M  = (sigma - delta*P)*(K+P)/(alpha*P);
        ev = eig(J(P,M));
        fprintf(' Coexist  E : P*=%.4f, M*=%.4f\n', P, M);
        fprintf('    eigenvalues = [% .4f , % .4f]\n', real(ev(1)), real(ev(2)));
        fprintf('    --> %s\n\n', classify(ev));
        Pstar_list(end+1) = P;
        Mstar_list(end+1) = M;
    end
end

f = @(t,y)[ sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1));
            theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2) ];

eps = 0.05;
fsp = @(t,y)[ (sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1)))/eps;
               theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2) ];

tspan = [0 200];  y0 = [5; 2];
[t1,Y1] = ode45 (f,   tspan, y0);
[t2,Y2] = ode15s(fsp, tspan, y0);

figure('Name','Time series','Color','w');
subplot(1,2,1);
plot(t1,Y1(:,1),'b','LineWidth',1.6); hold on;
plot(t1,Y1(:,2),'r','LineWidth',1.6);
xlabel('time t'); ylabel('concentration');
legend('Pollutant P(t)','Microbe M(t)','Location','best');
title('Full model (\epsilon = 1)'); grid on;

subplot(1,2,2);
plot(t2,Y2(:,1),'b','LineWidth',1.6); hold on;
plot(t2,Y2(:,2),'r','LineWidth',1.6);
xlabel('time t'); ylabel('concentration');
legend('Pollutant P(t)','Microbe M(t)','Location','best');
title('Singularly perturbed (\epsilon = 0.05)'); grid on;

figure('Name','Phase portrait','Color','w'); hold on;

for Pi = 0.5:1.5:11
    for Mi = 0.2:1.0:4.5
        [~,Y] = ode45(f,[0 400],[Pi;Mi]);
        plot(Y(:,1),Y(:,2),'Color',[0.7 0.7 0.7]);
    end
end

Pg = linspace(0.05,11,600);
Mnull = (sigma - delta*Pg).*(K+Pg)./(alpha*Pg);
plot(Pg, max(Mnull,0), 'b', 'LineWidth', 2);

for k = 1:numel(Pstar_list)
    xline(Pstar_list(k), 'g--', 'LineWidth', 2);
end

plot(P0,0,'ks','MarkerFaceColor','k','MarkerSize',9);
for k = 1:numel(Pstar_list)
    plot(Pstar_list(k),Mstar_list(k),'ro', ...
         'MarkerFaceColor','r','MarkerSize',9);
end

xlabel('Pollutant  P'); ylabel('Microbe  M');
title('Phase portrait: nullclines, equilibria, bistability');
legend('trajectories','dP/dt = 0 nullcline','dM/dt = 0 nullcline', ...
       'Location','northeast');
axis([0 11 0 5]); grid on; box on;

gam = linspace(0, 0.25, 800);
figure('Name','Bifurcation','Color','w'); hold on;
for g = gam
    c = [g, mu + g*K - theta, mu*K];
    r = roots(c);
    for i = 1:numel(r)
        P = r(i);
        if isreal(P) && P>0 && P < sigma/delta
            M = (sigma - delta*P)*(K+P)/(alpha*P);
            detJ = (alpha*P/(K+P))*M*(theta*K/(K+P)^2 - g);
            if detJ > 0
                plot(g, M, 'b.', 'MarkerSize', 6);
            else
                plot(g, M, 'r.', 'MarkerSize', 6);
            end
        end
    end
end
xlabel('Toxicity coefficient  \gamma');
ylabel('Equilibrium biomass  M^*');
title('Saddle-node bifurcation (blue = stable, red = saddle)');
grid on; box on;

fprintf('============================================================\n');
fprintf(' Done. Check Figures 1-3 and copy the console output above.\n');
fprintf('============================================================\n');

function s = classify(ev)
    re = real(ev);
    if all(re < 0)
        s = 'STABLE (sink)';
    elseif all(re > 0)
        s = 'UNSTABLE (source)';
    elseif any(re > 0) && any(re < 0)
        s = 'SADDLE';
    else
        s = 'NON-HYPERBOLIC (borderline)';
    end
end

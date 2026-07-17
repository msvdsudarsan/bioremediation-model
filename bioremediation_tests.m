clear; clc; close all;

sigma = 1.0; delta = 0.1; alpha = 0.8; K = 1.0;
theta = 0.6; mu = 0.2; gamma = 0.05;
opt = odeset('RelTol',1e-9,'AbsTol',1e-11);

disp('============================================================');
disp(' SUPPLEMENTARY TESTS');
disp('============================================================');

P0 = 5; M0 = 2; Tend = 200; T0 = 20;
tgrid = linspace(T0,Tend,2000);
Pslow = @(M) ((sigma-delta*K-alpha*M) + sqrt((sigma-delta*K-alpha*M).^2 + 4*delta*sigma*K))/(2*delta);
slow  = @(t,M) theta*Pslow(M)*M/(K+Pslow(M)) - mu*M - gamma*Pslow(M)*M;
[tr,Mr] = ode15s(slow,[0 Tend],M0,opt);
Mr_i = interp1(tr,Mr,tgrid,'pchip');
fullsys = @(t,y,ep) [ (1/ep)*(sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1))); ...
                      theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2) ];
epsv = [1 0.5 0.2 0.1 0.05 0.02 0.01];
err  = zeros(size(epsv));
for k = 1:numel(epsv)
    [tf,yf] = ode15s(@(t,y) fullsys(t,y,epsv(k)),[0 Tend],[P0;M0],opt);
    err(k) = max(abs(interp1(tf,yf(:,2),tgrid,'pchip') - Mr_i));
end
pc = polyfit(log(epsv),log(err),1);
fprintf('\nTEST A - SP convergence (t>=%g):\n',T0);
fprintf('   eps      max error\n');
for k = 1:numel(epsv)
    fprintf('  %5.3f   %10.3e\n',epsv(k),err(k));
end
fprintf('   fitted log-log slope = %.3f   (=> error ~ O(eps^%.2f))\n',pc(1),pc(1));

figure('Color','w');
loglog(epsv,err,'o-','LineWidth',1.5); grid on; hold on;
loglog(epsv,err(end)/epsv(end)*epsv,'k--');
xlabel('epsilon'); ylabel('max_t | M full - M reduced |');
title('Test A: SP convergence (slope ~ 1 => O(eps))');
legend('measured','slope-1 ref','Location','northwest');
try
    exportgraphics(gcf,'fig_sp_convergence.png','Resolution',150);
    disp('   [saved fig_sp_convergence.png]');
catch
    try, saveas(gcf,'fig_sp_convergence.png'); disp('   [saved via saveas]'); catch, disp('   [save failed - use Figure > Save As]'); end
end

Kv = linspace(0.5,2.0,200);
gc_closed = (sqrt(theta)-sqrt(mu)).^2 ./ Kv;
gc_num = zeros(size(Kv));
for i = 1:numel(Kv)
    gg = linspace(1e-4,0.3,60000);
    bb = mu + gg*Kv(i) - theta;
    disc = bb.^2 - 4*gg*(mu*Kv(i));
    gc_num(i) = gg(find(disc>=0,1,'last'));
end
fprintf('\nTEST B - fold boundary gamma_c(K):\n');
fprintf('   max |closed form - numerical| over K = %.3e\n',max(abs(gc_closed-gc_num)));
fprintf('   nominal gamma_c(K=1) = %.5f\n',(sqrt(theta)-sqrt(mu))^2/K);

figure('Color','w');
plot(Kv,gc_closed,'b-','LineWidth',2); hold on;
plot(Kv,gc_num,'r--','LineWidth',1.2);
plot(K,(sqrt(theta)-sqrt(mu))^2/K,'k*','MarkerSize',10);
xlabel('half-saturation K'); ylabel('critical toxicity gamma_c'); grid on;
legend('closed form','numerical fold','nominal K=1');
title('Test B: tipping boundary gamma_c(K)');
try
    exportgraphics(gcf,'fig_tipping_boundary.png','Resolution',150);
    disp('   [saved fig_tipping_boundary.png]');
catch
    try, saveas(gcf,'fig_tipping_boundary.png'); disp('   [saved via saveas]'); catch, disp('   [save failed - use Figure > Save As]'); end
end

E1 = [0.6277,3.0379]; E0 = [10,0];
Pg = linspace(0.2,10,40); Mg = linspace(0,5,35);
lab = zeros(numel(Mg),numel(Pg));
f1 = @(t,y) [ sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1)); ...
              theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2) ];
for a = 1:numel(Pg)
    for c = 1:numel(Mg)
        [~,y] = ode45(f1,[0 400],[Pg(a);Mg(c)]);
        e = y(end,:);
        lab(c,a) = (norm(e-E1) < norm(e-E0));
    end
end
frac = mean(lab(:));
fprintf('\nTEST C - basin of attraction:\n');
fprintf('   fraction of initial states -> recovery (E1) = %.3f\n',frac);

figure('Color','w');
imagesc(Pg,Mg,lab); set(gca,'YDir','normal');
colormap([0.85 0.3 0.3; 0.3 0.5 0.85]); hold on;
plot(E1(1),E1(2),'w.','MarkerSize',22);
plot(E0(1),E0(2),'ks','MarkerFaceColor','k','MarkerSize',9);
xlabel('initial pollutant P0'); ylabel('initial biomass M0');
title(sprintf('Test C: basins (blue=recovery, red=washout), recovery frac=%.2f',frac));
try
    exportgraphics(gcf,'fig_basins.png','Resolution',150);
    disp('   [saved fig_basins.png]');
catch
    try, saveas(gcf,'fig_basins.png'); disp('   [saved via saveas]'); catch, disp('   [save failed - use Figure > Save As]'); end
end

disp(' ');
disp('============================================================');
disp(' Done. Post the console block + attach the 3 PNG figures.');
disp('============================================================');

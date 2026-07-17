clear; clc;

sigma = 1.0; delta = 0.1; alpha = 0.8; K = 1.0;
theta = 0.6; mu = 0.2; gamma = 0.05;
KT    = 1.0;

fLin  = @(t,y) [ sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1));
                 theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2) ];
fHill = @(t,y) [ sigma - delta*y(1) - alpha*y(1)*y(2)/(K+y(1));
                 theta*y(1)*y(2)/(K+y(1)) - mu*y(2) - gamma*y(1)*y(2)/(KT+y(1)) ];

IC = [ 5.0 2.0;
       9.0 0.2;
       2.0 3.0;
       1.0 4.0;
       8.0 0.5 ];

fprintf('model   start      P*        M*\n');
for k = 1:size(IC,1)
    [~,yL] = ode45(fLin,  [0 600], IC(k,:));
    [~,yH] = ode45(fHill, [0 600], IC(k,:));
    fprintf('linear  %2d   %8.4f %8.4f\n', k, yL(end,1), yL(end,2));
    fprintf('hill    %2d   %8.4f %8.4f\n', k, yH(end,1), yH(end,2));
end

fracLin  = basin_fraction(fLin,  [0.6277 3.0379]);
fracHill = basin_fraction(fHill, [0.5715 3.2410]);
fprintf('\nrecovery basin fraction (linear): %.3f\n', fracLin);
fprintf('recovery basin fraction (hill)  : %.3f\n', fracHill);

P0      = sigma/delta;
lamLin  = theta*P0/(K+P0) - mu - gamma*P0;
lamHill = theta*P0/(K+P0) - mu - gamma*P0/(KT+P0);
fprintf('\nwashout E0 transverse eigenvalue (linear): %+.4f\n', lamLin);
fprintf('washout E0 transverse eigenvalue (hill)  : %+.4f\n', lamHill);

function frac = basin_fraction(f, Erec)
    E0 = [10 0];
    Pg = linspace(0.2, 10, 40);
    Mg = linspace(0,   5, 35);
    lab = zeros(numel(Mg), numel(Pg));
    for a = 1:numel(Pg)
        for c = 1:numel(Mg)
            [~,y] = ode45(f, [0 400], [Pg(a); Mg(c)]);
            e = y(end,:);
            lab(c,a) = (norm(e - Erec) < norm(e - E0));
        end
    end
    frac = mean(lab(:));
end

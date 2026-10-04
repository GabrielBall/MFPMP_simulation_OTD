function plot_solution(cfg,Y,Rho,u)

% Density
figure('Name','Mean-field density');
imagesc(cfg.time,cfg.cellCentres,Rho);
axis xy tight;
hold on;
plot(cfg.time,Y,'w','LineWidth',2);
yline(cfg.target,'w--','LineWidth',1.5);
xlabel('t');
ylabel('x');
colorbar;

% End states
figure('Name','Initial and final density');
plot(cfg.cellCentres,Rho(:,1),'LineWidth',1.5);
hold on;
plot(cfg.cellCentres,Rho(:,end),'LineWidth',1.5);
xline(cfg.target,'--','LineWidth',1.5);
xlabel('x');
ylabel('\rho');
legend('Initial','Final','Target','Location','best');
grid on;

% Control
figure('Name','Prescribed control');
stairs(cfg.time(1:end-1),u,'LineWidth',1.5);
xlabel('t');
ylabel('u(t)');
grid on;

end

function plot_optimisation_results(cfg,Y,Rho,u,costHistory,residualHistory)

if nargin < 5
    costHistory = [];
end

if nargin < 6
    residualHistory = [];
end

% Density
figure('Name','Optimised mean-field density');
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
figure('Name','Optimised control');
stairs(cfg.time(1:end-1),u,'LineWidth',1.5);
xlabel('t');
ylabel('u(t)');
grid on;

% Histories
if ~isempty(costHistory)
    figure('Name','Objective history');
    semilogy(0:numel(costHistory)-1,costHistory,'o-');
    xlabel('Iteration');
    ylabel('J');
    grid on;
end

if ~isempty(residualHistory)
    figure('Name','PMP residual history');
    semilogy(1:numel(residualHistory),residualHistory,'o-');
    xlabel('Iteration');
    ylabel('PMP residual');
    grid on;
end

end

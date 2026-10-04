function [value,gradRho,gradY] = binary_voter(~,rho,y,cfg,parameters)

rho = rho(:);

threshold = parameters.threshold;
transitionWidth = parameters.transitionWidth;
weight = parameters.weight;

vote = 0.5.*(1+tanh((cfg.cellCentres(:)-threshold)./transitionWidth));

value = weight.*cfg.dx.*sum(vote.*rho);
gradRho = weight.*cfg.dx.*vote;
gradY = zeros(size(y));

end
function [leftCoefficient,rightCoefficient,dleftDv,drightDv] = mean_field_flux_coefficients(velocity,diffusion,dx)

velocity = velocity(:);

% Zero diffusion
if diffusion == 0
    leftCoefficient = max(velocity,0);
    rightCoefficient = max(-velocity,0);

    tolerance = 100.*eps(max(1,max(abs(velocity))));
    positive = velocity > tolerance;
    negative = velocity < -tolerance;
    zeroVelocity = ~(positive | negative);

    dleftDv = zeros(size(velocity));
    drightDv = zeros(size(velocity));
    dleftDv(positive) = 1;
    drightDv(negative) = -1;
    dleftDv(zeroVelocity) = 0.5;
    drightDv(zeroVelocity) = -0.5;
    return;
end

% Chang-Cooper coefficients
pecletNumber = velocity.*dx./diffusion;
leftCoefficient = (diffusion./dx).*bernoulli_function(-pecletNumber,false);
rightCoefficient = (diffusion./dx).*bernoulli_function(pecletNumber,false);
dleftDv = -bernoulli_function(-pecletNumber,true);
drightDv = bernoulli_function(pecletNumber,true);

end

function value = bernoulli_function(z,derivative)

value = zeros(size(z));
small = abs(z) < 1e-6;
largePositive = z > 50;
largeNegative = z < -50;
regular = ~(small | largePositive | largeNegative);

% Small argument
if any(small)
    zs = z(small);
    if derivative
        value(small) = -0.5+(1/6).*zs-(1/180).*zs.^3+(1/5040).*zs.^5;
    else
        value(small) = 1-0.5.*zs+(1/12).*zs.^2-(1/720).*zs.^4;
    end
end

% Regular argument
if any(regular)
    zr = z(regular);
    denominator = expm1(zr);
    if derivative
        value(regular) = (denominator-zr.*(denominator+1))./(denominator.^2);
    else
        value(regular) = zr./denominator;
    end
end

% Large positive
if any(largePositive)
    zp = z(largePositive);
    if derivative
        value(largePositive) = (1-zp).*exp(-zp);
    else
        value(largePositive) = zp.*exp(-zp);
    end
end

% Large negative
if any(largeNegative)
    if derivative
        value(largeNegative) = -1;
    else
        value(largeNegative) = -z(largeNegative);
    end
end

end

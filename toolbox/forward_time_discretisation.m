function zNext = forward_time_discretisation(z,h,scheme,A,f,args,argsNext)

z = z(:);
N = numel(z);

if isempty(args)
    args = {};
end
if isempty(argsNext)
    argsNext = args;
end

% Current stage
A = A(z,args{:});
if isempty(f)
    f = zeros(N,1);
else
    f = f(z,args{:});
end

% Scheme
switch lower(char(scheme))
    case 'explicit'
        zNext = z+h.*(A*z+f);

    case 'semi_implicit'
        zNext = (speye(N)-h.*A)\(z+h.*f);

    case 'imex'
        zTilde = forward_time_discretisation(z,h,'explicit',A,f,args,args);
        ATilde = A(zTilde,argsNext{:});

        if isempty(f)
            fTilde = zeros(N,1);
        else
            fTilde = f(zTilde,argsNext{:});
        end

        zNext = (speye(N)-0.5.*h.*ATilde)\(z+0.5.*h.*(A*z+f+fTilde));

    otherwise
        error('Unknown time discretisation scheme: %s.',scheme);
end

end

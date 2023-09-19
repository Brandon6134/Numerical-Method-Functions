function [du] = richardson21( D, u, x, h, N_max, eps_abs )
% Argument Checking
if ~ isa(D,'function_handle')
    throw(MException('MATLAB:invalid_argument', 'The argument D is not a valid function'));

end

if ~ isa(u,'function_handle')
    throw(MException('MATLAB:invalid_argument', 'The argument u is not a valid function'));

end

if ~ isscalar(x)
    throw(MException('MATLAB:invalid_argument', 'The argument x is not scalar'));

end

if ~ isscalar(h)
    throw(MException('MATLAB:invalid_argument', 'The argument h is not scalar'));

end

if ~ isscalar(N_max) || (N_max ~= round(N_max))
    throw(MException('MATLAB:invalid_argument', 'The argument N_max is not positive integer'));

end

if ~ isscalar(eps_abs) || (eps_abs <= 0)
    throw(MException('MATLAB:invalid_argument', 'The argument eps_abs is not positive integer'));

end

% creates the initial square matrix and sets element value to R(1,1)
R = zeros(N_max+1,N_max+1);
R(1,1) = D(u,x,h);

%creates the loops that implements the richardson extrapolation method, and returns a final value
for i = 1: N_max
    R(i+1,1) = D(u,x,h/(2^i));
    for j = 1 : i
        R(i+1,j+1)=(2^(j+1)*R(i+1,j)-R(i,j))/(2^(j+1)-1);
    end

    if (abs(R(i+1,i+1)-R(i,i))<eps_abs)
        du = R(i+1,i+1);
        break;
    elseif(i == N_max)
        error('MATLAB:invalid_argument', 'Richardson extrapolation did not converge');
    end
end
end


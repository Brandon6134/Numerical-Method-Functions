% Secant method
%
% The secant function uses the secant method to find the root of the given
% function. Using the 2 initial x values, it searches for a root given the
% parameters of how accurate the user wants the root to be.

% Getting into specifics, it solves for the inital root based on the
% parameters x1 and x0, then checks if the difference between this root and
% the previous root is less than eps_step, and if its close enough to 0
% according eps_abs. If this condition is not met, x0 and x1 are given new
% values based on the previous root search then the function loops until
% the 2 conditions with eps_step and eps_abs are true. If it cannot find a
% suitable root within the given amount of loops desired (N_max), function
% gives an error.
%
% Parameters                 
% ==========                 
%    f         The function that the secant function tries to find the root of.
%    x0        Initial leftmost x value boundary, controls where root search begins
%    x1        Initial rightmost x value boundary, controls where root search begins
%    eps_step  Roots found must have differences lower this variable; this is one condition that stops root search loop.
%    eps_abs   Root found must have a y value lower than this variable; this is one condition that stops root search loop.
%    N_max     Maximum number of loops that the root search can undergo.
%
% Return Values              
% =============
%    x2        Approximation of the root found for the given function.

function [x2] = secant( f, x0, x1, eps_step, eps_abs, N_max )
    for i = 1:N_max
        x2 = (f(x0)*(x1-x0)*x1 + f(x1)*(x0-x1)*x0)/( (f(x0)*(x1-x0)+ f(x1)*(x0-x1)));
        %above line solves for root
        
        %below if statement checks if root differences are below desired
        %value and if the root found gives a value close enough to zero
        %based on the given eps_abs value.
        if abs(x2-x1) < eps_step && abs(f(x2)) < eps_abs
            return;
        end

        x0 = x1;
        x1 = x2;
    end
    throw(MException('MATLAB:numeric_exception','convergence did not occur'));
    %gives error if exceeds N_max number of loops
end


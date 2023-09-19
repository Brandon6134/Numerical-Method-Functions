function [t_out, y_out] = dp45( f, t_rng, y0, h, eps_abs )
% Dormand Prince
%
% Using Dormand-Prince’s algorithm to minimise errors and improve
% accuracy of the approximations of 1st and 2nd order initial-value
% problems.
%
% Parameters
% ==========
% f : a function handle to the bivariate function f(t,y). Provides
% slope
% t_rng : a row vector of two values [t0, tfinal]
% y0 : the initial condition
% h : a very small change
% eps_abs: absolute error
%
% Return Values
% =============
% t_out : a row vector of n equally spaced values from t0
% to tfinal
% y_out : a row vector of n values where the first
% point(y_out(1)) equals y0 and the k
% point approximates y(t) at the t_out(k) from 2 to n
%
% Argument Checking
if ~isa( f, 'function_handle' )
throw( MException( 'MATLAB:invalid_argument', ...
'the argument f is not a function handle' ) );
end
if ~all( size( t_rng ) == [1, 2] )
throw( MException( 'MATLAB:invalid_argument', ...
'the argument t_rng is not a row vector with two entries' ) );
end
if ~isscalar( y0 )
throw( MException( 'MATLAB:invalid_argument', ...
'the argument y0 is not a scalar' ) );
end
if ~isscalar(h)
throw( MException( 'MATLAB:invalid_argument',...
'the argument h is not a positive scalar' ) );
end
if ~isscalar(eps_abs) || (eps_abs <= 0)
throw( MException( 'MATLAB:invalid_argument',...
'the argument eps_abs is not a positive integer' ) );
end
% Initialize t_out and y_out
A = [0 0 0 0 0 0 0;
1 0 0 0 0 0 0;
1/4 3/4 0 0 0 0 0;
11/9 -14/3 40/9 0 0 0 0;
4843/1458 -3170/243 8056/729 -53/162 0 0 0;
9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
35/384 0 500/1113 125/192 -2187/6784 11/84 0]';

by = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]';
bz = [35/384 0 500/1113 125/192 -2187/6784 11/84 0]';
c= [0 1/5 3/10 4/5 8/9 1 1]';

% Initialize the location to k = 1
k = 1;
t0 = t_rng(1);
tf = t_rng(2);
t_out(k) = t0;
y_out(k) = y0;
% Looping
% Use Dormand Prince to find two approximations,
% y_tmp and z_tmp to approximate y(t) at
% t = t_out(k) + h for the current value of h
while t_out(k) < tf
K = zeros(1,7);
for m = 1:7
K(m) = f(t_out(k) + h*c(m) , y_out(k) + h*c(m)*K*A(:,m));
end

y_tmp = y_out(k) + h*K*by;
z_tmp = y_out(k) + h*K*bz;

% Calculate the scaling factor ‘s’
s = ((eps_abs*h)/(2*(tf-t0)*abs(y_tmp - z_tmp)))^(1/4);


% if s >= 2,
% We use z_tmp to approximate y_out(k + 1)
% t_out(k + 1) is the previous t-value plus h
% Increment k and double the value of h for the
% next iteration
if s >= 2
y_out(k+1) = z_tmp;
t_out(k+1) = t_out(k) + h;
k = k+1;
h = h*2;

% else if s >= 1,
% We use z_tmp to approximate y_out(k + 1)
% t_out(k + 1) is the previous t-value plus h
% In this case, h is neither too large or too
% small, so only increment k
elseif s >= 1
y_out(k+1) = z_tmp;
t_out(k+1) = t_out(k) + h;
k = k+1;


% else s < 1
% Divide h by two and try again with the smaller
% value of h (just go through the loop again
% without updating t_out, y_out, or k)
elseif s < 1
    h= h/2;
end

% We must make one final check before we end the loop:
% if t_out(k) + h > tf, we must reduce the
% size of h so that t_out(k) + h == tf
if (t_out(k) + h) > tf

h = tf - t_out(k);
end
end
end


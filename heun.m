function [ t_out, y_out ] = heun( f, t_rng, y0, n )
t0 = t_rng(1,1);
t_final = t_rng(1,2);
h = (t_final - t0)/(n-1);
t_out = linspace(t0, t_final, n);
y_out(1) = y0;
 
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
if ~isscalar( n ) || (n <= 0) || (n ~= round( n ))
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n is not a positive integer' ) );
end
 
for k = 1:(n-1)
    K1 = f(t_out(k),y_out(k));
    K2 = f(t_out(k) + h, y_out(k)+ h.*K1);
    y_out(k+1) = y_out(k)+h.*((K1+K2)/2);
end



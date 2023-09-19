function [x_out, u_out] = shooting( s1, s2, f, x_rng, u_bndry, h, eps_abs, eps_step, N_max )
% Argument Checking
if ~isa( f, 'function_handle' )
throw( MException( 'MATLAB:invalid_argument', ...
'the argument f is not a function handle' ) );
end
if ~all( size( t_rng ) == [1, 2] )
throw( MException( 'MATLAB:invalid_argument', ...
'the argument t_rng is not a row vector with two entries' ) );
end
if ~isscalar(h)
throw( MException( 'MATLAB:invalid_argument',...
'the argument h is not a positive scalar' ) );
end
if ~isscalar(eps_abs) || (eps_abs <= 0)
throw( MException( 'MATLAB:invalid_argument',...
'the argument eps_abs is not a positive integer' ) );
end

end


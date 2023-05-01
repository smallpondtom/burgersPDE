function burgers_solution_test ( )

%*****************************************************************************80
%
%% burgers_solution_test() tests burgers_solution().
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 December 2018
%
%  Author:
%
%    John Burkardt
%
  addpath ( '../burgers_solution' )

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'burgers_solution_test():\n' );
  fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Test burgers_solution().\n' );

  burgers_solution_test01 ( );
  pause ( 5 );
  burgers_solution_test02 ( );
  pause ( 5 );
  burgers_solution_test03 ( );
  pause ( 5 );
  burgers_solution_test04 ( );
  pause ( 5 );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'burgers_solution_test():\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../burgers_solution' )

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% timestamp() prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end


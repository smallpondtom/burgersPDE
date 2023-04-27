function burgers_solution_test03 ( )

%*****************************************************************************80
%
%% burgers_solution_test03() sets up a small test case for solution #2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 September 2015
%
%  Author:
%
%    John Burkardt
%
  tn = 11;
  xn = 11;
  nu = 0.5;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'burgers_solution_test03()\n' );
  fprintf ( 1, '  Compute analytic solution #2 to the Burgers equation.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Viscosity NU = %g\n', nu );
  fprintf ( 1, '  NX = %d\n', xn );
  fprintf ( 1, '  NT = %d\n', tn );

  xlo = 0.0;
  xhi = 2.0 * pi;
  x = linspace ( xlo, xhi, xn );
  r8vec_print ( xn, x, '  X grid points:' );

  tlo = 0.0;
  thi = 1.0;
  t = linspace ( tlo, thi, tn );
  r8vec_print ( tn, t, '  T grid points:' );

  u = burgers_viscous_time_exact2 ( nu, xn, x, tn, t );

  r8mat_print ( xn, tn, u, '  U(X,T) at grid points:' );

  surf ( x, t, u' );
  xlabel ( '<--- X --->' );
  ylabel ( '<--- T --->' );
  zlabel ( '<--- U(X,T) --->' );
  title ( 'burgers\_solution\_test03' )

  filename = 'burgers_solution_test03.txt';

  r8mat_write ( filename, xn, tn, u );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Data written to file "%s"\n', filename );

  return
end

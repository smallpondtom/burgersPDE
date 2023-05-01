function burgers_solution_test01 ( )

%*****************************************************************************80
%
%% burgers_solution_test01() tests sets up a small test case.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    17 November 2011
%
%  Author:
%
%    John Burkardt
%
  vtn = 11;
  vxn = 11;
  nu = 0.01 / pi;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'burgers_solution_test01():\n' );
  fprintf ( 1, '  Compute an analytic solution to the Burgers equation.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Viscosity NU = %g\n', nu );
  fprintf ( 1, '  NX = %d\n', vxn );
  fprintf ( 1, '  NT = %d\n', vtn );

  xlo = -1.0;
  xhi = +1.0;
  vx = linspace ( xlo, xhi, vxn );
  r8vec_print ( vxn, vx, '  X grid points:' );

  tlo = 0.0;
  thi = 3.0 / pi;
  vt = linspace ( tlo, thi, vtn );
  r8vec_print ( vtn, vt, '  T grid points:' );

  vu = burgers_viscous_time_exact1 ( nu, vxn, vx, vtn, vt );

  r8mat_print ( vxn, vtn, vu, '  U(X,T) at grid points:' );

  surf ( vx, vt, vu' );
  xlabel ( '<--- X --->' );
  ylabel ( '<--- T --->' );
  zlabel ( '<--- U(X,T) --->' );
  title ( 'burgers\_solution\_test01' );

  filename = 'burgers_solution_test01.txt';

  r8mat_write ( filename, vxn, vtn, vu );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Data written to file "%s"\n', filename );

  return
end

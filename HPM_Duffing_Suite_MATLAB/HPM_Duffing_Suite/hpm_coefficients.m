function [U, XF, wF] = hpm_coefficients()
%HPM_COEFFICIENTS  Exact rational coefficients of the HPM solutions.
%
%  Non-dimensional problem:  X'' + X + eps*X^3 = 0,  X(0)=1, X'(0)=0,
%  with tau = sqrt(C)*t and eps = D*A^2/C (A = initial amplitude, here 1).
%
%  U  : plain HPM,  X = sum_n eps^n u_n(tau).  Each row of U is
%       [n  numerator  denominator  k  m  type]  meaning
%       (num/den) * tau^k * cos(m*tau)   (type = 0)
%       (num/den) * tau^k * sin(m*tau)   (type = 1)
%  XF : HPM with frequency expansion (HPM-FE), stretched time s = Omega*tau,
%       X = sum_n eps^n X_n(s); row [n num den m] means (num/den)*cos(m*s).
%  wF : Omega^2 = 1 + sum_j wF(j)*eps^j  (j = 1..4).
%
%  The coefficients were generated symbolically and are verified by
%  verify_terms.m against direct numerical integration of the recurrences.
U = [
  0 1 1 0 1 0;
  1 -1 32 0 1 0;
  1 -3 8 1 1 1;
  1 1 32 0 3 0;
  2 23 1024 0 1 0;
  2 -9 128 2 1 0;
  2 3 32 1 1 1;
  2 -3 128 0 3 0;
  2 -9 256 1 3 1;
  2 1 1024 0 5 0;
  3 -547 32768 0 1 0;
  3 135 4096 2 1 0;
  3 -207 4096 1 1 1;
  3 9 1024 3 1 1;
  3 297 16384 0 3 0;
  3 -81 4096 2 3 0;
  3 279 8192 1 3 1;
  3 -3 2048 0 5 0;
  3 -15 8192 1 5 1;
  3 1 32768 0 7 0;
  4 6713 524288 0 1 0;
  4 -1359 65536 2 1 0;
  4 27 32768 4 1 0;
  4 8997 262144 1 1 1;
  4 -99 16384 3 1 1;
  4 -15121 1048576 0 3 0;
  4 1539 65536 2 3 0;
  4 -3915 131072 1 3 1;
  4 243 32768 3 3 1;
  4 883 524288 0 5 0;
  4 -225 131072 2 5 0;
  4 825 262144 1 5 1;
  4 -9 131072 0 7 0;
  4 -21 262144 1 7 1;
  4 1 1048576 0 9 0;
];
XF = [
  0 1 1 1;
  1 -1 32 1;
  1 1 32 3;
  2 23 1024 1;
  2 -3 128 3;
  2 1 1024 5;
  3 -547 32768 1;
  3 297 16384 3;
  3 -3 2048 5;
  3 1 32768 7;
  4 6713 524288 1;
  4 -15121 1048576 3;
  4 883 524288 5;
  4 -9 131072 7;
  4 1 1048576 9;
];
wF = [3/4, -3/128, 9/512, -1779/131072];
end

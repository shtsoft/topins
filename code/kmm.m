%{
% kmm stands for kane-mele model
%
% kmm.m is a matlab/octave script that
% - plots the energy bands of the model on a fundamental domain
% - calculates the chern numbers C of the model
% - calculates the thermodynamic pressures p and K of the model and confirm the conjecture p = K
%
% the hilbert space of the kane-mele model can be identified with $\mathbb{C}^{16m^{2}}$ if we 
% consider a finite honeycomb crystal with periodic boundary conditions
%
% kmm takes as arguments (in the above context):
% - m: m^(2) is the number of crystal snippets  
%      0 B A 0
%      A 0 0 B
%      0 B A 0
%      A 0 0 B
% - d: dist(A, B) - the nearest neighbor distance
% - n: (n + 1) * (n + 1) is the number of grid points in the fundamental domain
% - tnn, tse, tso and tra: prefactors of the corresponding hamiltonian terms
% - Bvar: the real valued variation parameter of an external magnetic field
% - Bint: the integer valued strength paramter of an external magnetic field
% - betakb: the inverse of temperature T
% - mu: the chemical potential
%
% kmm returns the positive difference of p and K
%
% WARNING:
% computation time and memory usage depends mainly on m, n where m = 8, ..., 12 and 
% n = 50, ..., 120 works pretty well
%}
function pdiffofpK = kmm(m, n, d, tnn, tse, tso, tra, Bvar, Bint, betakb, mu)
  a = sqrt(3) * d;    % lattice constant
  b = 4 * pi / (a * sqrt(3));    % reciprocal lattice constant


  B = (2 * pi / (a * a * abs(sin(2 * pi / 3)))) * Bint + Bvar;    % ext mag. field B


  initdat(m, d, tnn, tse, tso, tra, B);


  hx = b / n;
  hy = b * cos(pi / 6) / n;
  kx = 0:hx:b;
  ky = 0:hy:b * cos(pi / 6);
  [Kx, Ky] = meshgrid(kx, ky);

  for j = 1:n + 1
    rx(j) = hx;
    ry(j) = hy;
  end

  rx(1) = hx / 2;
  ry(1) = hy / 2;
  rx(n + 1) = hx / 2;
  ry(n + 1) = hy / 2;


  dimhf = 4;    % = rows(Hf(0, 0)) in octave;
  C(1:dimhf) = 0;
  K(1:dimhf) = 0;

  for jx = 1:n + 1
    for jy = 1:n + 1
      [V, D] = eig(Hf(kx(jx), ky(jy)));
      [Vhx1, Dhx1] = eig(Hf(kx(jx) + hx, ky(jy)));
      [Vhx2, Dhx2] = eig(Hf(kx(jx) - hx, ky(jy)));
      [Vhy1, Dhy1] = eig(Hf(kx(jx), ky(jy) + hy));
      [Vhy2, Dhy2] = eig(Hf(kx(jx), ky(jy) - hy));
      for l = 1:dimhf
        P = V(:, l) * V(:, l)';
        DxP = (Vhx1(:, l) * Vhx1(:, l)' - Vhx2(:, l) * Vhx2(:, l)') / (2 * hx);
        DyP = (Vhy1(:, l) * Vhy1(:, l)' - Vhy2(:, l) * Vhy2(:, l)') / (2 * hy);
        TrP = trace(P * DxP * DyP - P * DyP * DxP);
        TrHf = trace(P * DxP * (Hf(kx(jx), ky(jy)) - D(l, l)) * DyP);
        Arg =  - betakb * (D(l, l) + Bvar * imag(TrHf) - mu);
        if Arg > 30
          K(l) = K(l)  +  rx(jx)  *  ry(jy)  *  (1  -  Bvar  *  1i  * TrP)  *  Arg;
        else
          K(l) = K(l) + rx(jx) * ry(jy) * (1 - Bvar * 1i * TrP) * log(1 + exp(Arg));
        end
        C(l) = C(l) + rx(jx) * ry(jy) * TrP;
        Z(jy, jx, l) = D(l, l);
      end
    end
  end

  hold off;
  for l = 1:dimhf
    S(:, :) = Z(:, :, l);
    surf(Kx, Ky, real(S));
    hold on;
  end


  global dimhx;
  global Hx;
  p = 0;
  Dmu = -betakb * eig(Hx - mu);
  for l = 1:dimhx
    if Dmu(l) > 30
      p = p + Dmu(l);
    else
      p = p + log(1 + exp(Dmu(l)));
    end
  end


  C = C / (2 * pi * 1i);
  p = p / (16 * m * m * betakb);
  K = K / (cos(pi / 6) * b * b * betakb);

  pdiffofpK = abs(p - (K(1) + K(2) + K(3) + K(4)) / 4);


  disp('chern numbers C(i):');
  disp(C);
  disp(' ');
  disp('thermopressure p:');
  disp(p);
  disp(' ');
  disp('thermopressure K(i):');
  disp(K);
  disp(' ');
  disp('|p - K| = |p - (sum_{i}K(i) / 4)|:');
  disp(pdiffofpK);
end



%{
% initialize hamiltonian prefactors globally
%
% initialize the lattice vectors a1, a2, a3 globally
%
% initialize the pauli matrices sx, sy, sz locally
% initialize some tensorproducts sxid, sxsx, ... globally
%
% initialize hilbert space dimension dimhx globally
% initialize nearest neighbor translation operators Tpd1, Tmd1, Tpd2, ... locally
% initialize next nearest neighbor translation ops TApa1, TAma1, TBpa1, ... locally
% initialize kane-mele hamiltonian nearest neighbor term hxnn locally
% initialize kane-mele hamiltonian self energy term hxse locally
% initialize kane-mele hamiltonian spin orbit term hxso locally
% initialize kane-mele hamiltonian rashba term hxra locally
% initialize kane-mele hamiltonian Hx globally
%}
function initdat(m, d, tnn, tse, tso, tra, B)
  global lamnn;
  global lamse;
  global lamso;
  global lamra;

  lamnn = tnn;
  lamse = tse;
  lamso = tso;
  lamra = tra;


  global a1x;
  global a1y;
  global a2x;
  global a2y;
  global a3x;
  global a3y;

  a1x = d * (3 / 2);
  a1y = d * (sqrt(3) / 2);
  a2x = -d * (3 / 2);
  a2y = d * (sqrt(3) / 2);
  a3x = 0;
  a3y = -d * sqrt(3);


  id = eye(2);
  sx = [0, 1; 1, 0];
  sy = [0, -1i; 1i, 0];
  sz = [1, 0; 0, -1];

  global sxid;
  global syid;
  global szid;
  global sxsx;
  global sysy;
  global szsz;
  global sxsy;
  global sysx;

  sxid = kron(sx, id);
  syid = kron(sy, id);
  szid = kron(sz, id);
  sxsx = kron(sx, sx);
  sysy = kron(sy, sy);
  szsz = kron(sz, sz);
  sxsy = kron(sx, sy);
  sysx = kron(sy, sx);


  global dimhx;
  dimhx = 16 * m * m;
  dimhxsl2 = dimhx / 2;
  dimhxsl4 = dimhx / 4;

  Tpd1 = zeros(dimhxsl2);
  Tpd2 = zeros(dimhxsl2);
  Tpd3 = zeros(dimhxsl2);
  for jpar = 1:4 * m
    x2 = (4 * m - jpar) * cos(pi / 6) * d;
    joff = 2 * m * (jpar - 1)
    if rem(jpar, 2) == 1
      for j = (2 + joff):2:(2 * m + joff)
        Tpd1(mod1(j + 2 * m, dimhxsl2), j) = exp(1i * (d * B * x2 / 2 - sqrt(3) * d * d * B / 8));
        Tpd2(mod1(j - 2 * m, dimhxsl2), j) = exp(1i * (d * B * x2 / 2 + sqrt(3) * d * d * B / 8));
        Tpd3(j - 1, j) = exp( - 1i * d * B * x2);
      end
    else
      for j = (1 + joff):2:(2 * m - 1 + joff)
        Tpd1(mod1(j + 2 * m, dimhxsl2), j) = exp(1i * (d * B * x2 / 2 - sqrt(3) * d * d * B / 8));
        Tpd2(mod1(j - 2 * m, dimhxsl2), j) = exp(1i * (d * B * x2 / 2 + sqrt(3) * d * d * B / 8));
        if j == 1 + joff
          Tpd3(j + 2 * m - 1, j) = exp( - 1i * d * B * x2);
        else
          Tpd3(j - 1, j) = exp( - 1i * d * B * x2);
        end
      end
    end
  end
  Tmd1 = conj(transpose(Tpd1));
  Tmd2 = conj(transpose(Tpd2));
  Tmd3 = conj(transpose(Tpd3));

  TApa1 = Tmd3 * Tpd2;
  TAma1 = Tmd2 * Tpd3;
  TBpa1 = Tpd2 * Tmd3;
  TBma1 = Tpd3 * Tmd2;
  TApa2 = Tmd1 * Tpd3;
  TAma2 = Tmd3 * Tpd1;
  TBpa2 = Tpd3 * Tmd1;
  TBma2 = Tpd1 * Tmd3;
  TApa3 = Tmd2 * Tpd1;
  TAma3 = Tmd1 * Tpd2;
  TBpa3 = Tpd1 * Tmd2;
  TBma3 = Tpd2 * Tmd1;

  hxnn = lamnn * ( ...
       kron(Tpd1 + Tmd1, id) ...
       + kron(Tpd2 + Tmd2, id) ...
       + kron(Tpd3 + Tmd3, id));

  hxse = eye(dimhxsl2);
  for j = 1:dimhxsl4
    hxse(2 * j - 1, 2 * j - 1) = -1;
  end
  hxse = lamse * kron(hxse, id);

  hxso = 1i * lamso * ( ...
       kron(TApa1 - TAma1, sz) ...
       - kron(TBpa1 - TBma1, sz) ...
       + kron(TApa2 - TAma2, sz) ...
       - kron(TBpa2 - TBma2, sz) ...
       + kron(TApa3 - TAma3, sz) ...
       - kron(TBpa3 - TBma3, sz));

  hxra = 1i * lamra * ( ...
       kron(Tpd1 - Tmd1,  -(sqrt(3) * sx + sy) / 2) ...
       + kron(Tpd2 - Tmd2, (sqrt(3) * sx - sy) / 2) ...
       + kron(Tpd3 - Tmd3, sy));

  global Hx;
  Hx = hxnn + hxse + hxso + hxra;
end



%{
% if f is the bloch-floquet transform then we have Hf(kx, ky) = f(Hx(f^(-1)(kx, ky)))
%}
function H = Hf(kx, ky)
  global lamnn;
  global lamse;
  global lamso;
  global lamra;
  global a1x;
  global a1y;
  global a2x;
  global a2y;
  global a3x;
  global a3y;
  global sxid;
  global syid;
  global szid;
  global sxsx;
  global sysy;
  global szsz;
  global sxsy;
  global sysx;

  cnn1k = 1 ...
        + cos(kx * a1x + ky * a1y) ...
        + cos(kx * a3x + ky * a3y);
  cnn2k = sin(kx * a1x + ky * a1y) ...
        - sin(kx * a3x + ky * a3y);
  cso1k = sin(kx * a1x + ky * a1y);
  cso2k = sin(kx * a2x + ky * a2y);
  cso3k = sin(kx * a3x + ky * a3y);
  cra1k = (sqrt(3) / 2) * (1 - cos(kx * a3x + ky * a3y));
  cra2k = cos(kx * a1x + ky * a1y) ...
        - (1 / 2) * (1 + cos(kx * a3x + ky * a3y));
  cra3k =  -(sqrt(3) / 2) * sin(kx * a3x + ky * a3y);
  cra4k =  -sin(kx * a1x + ky * a1y) ...
        - (1 / 2) * (sin(kx * a3x + ky * a3y));

  hfnn = lamnn * (cnn1k * sxid + cnn2k * syid);
  hfse = lamse * szid;
  hfso =  -2 * lamso * (cso1k + cso2k + cso3k) * szsz;
  hfra = lamra * (cra1k * sysx + cra2k * sysy + cra3k * sxsx + cra4k * sxsy)

  H = hfnn + hfse + hfso + hfra;
end



%{
% internal mod(a, m) is not suited to index a matrix so we use a slightly different version 
% mod1(a, m) which returns $b = m$ instead of $b = 0$ for the zero class 
%}
function b = mod1(a, m)
  b = a - m. * floor(a. / m)
  if b =  = 0
    b = m;
  end
end



%{
% sortinc sorts the solution of [V, D] = eig(H) in increasing order this means the eigenvalues of H
% will be stored increasing in D and the columns in V will be rearranged with respect to D
%}
function [V, D] = sortinc(V, D)
  [d, perm] = sort(diag(D));
  D = diag(d);
  V = V(:, perm)
end

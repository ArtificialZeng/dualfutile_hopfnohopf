clc;

syms v6 lambda real
% ============================================================
% Define G_red
% ============================================================
Gred = [ ...
0,          0,     -1,          0,          1,  0;
0,          0,      0,          v6,          0, -1;
1,          0,     -4,         -1,          0,  0;
-1,         -1, -40000, -(40002+v6),         -1,  0;
-sym(1)/100, -sym(1)/100, 0, -sym(1)/100, -sym(301)/100, -1;
0,          1,      0,          0,         -1, -4 ...
];
% ============================================================
% 1. Characteristic polynomial
% ============================================================
p = expand(det(lambda*eye(6) - Gred));
fprintf('====================================================\n');
fprintf('CHARACTERISTIC POLYNOMIAL\n');
fprintf('====================================================\n\n');
disp(collect(p,lambda))
% Constant coefficient
p0 = simplify(subs(p,lambda,0));
fprintf('\ndet Gred(v6) =\n');
disp(p0)
fprintf('\nValue of v6 for which det Gred(v6)=0:\n');
v6zero = solve(p0 == 0,v6);
disp(vpa(v6zero,30))
% ============================================================
% 2. Evaluate at v6 = 5 and v6 = 40000
% ============================================================
v6vals = [5, 40000];
tol = 1e-10;
for k = 1:length(v6vals)
v6v = v6vals(k);
fprintf('\n\n====================================================\n');
fprintf('v6 = %g\n',v6v);
fprintf('====================================================\n');
% --------------------------------------------------------
% Exact symbolic rank
% --------------------------------------------------------
Asym = subs(Gred,v6,v6v);
fprintf('Rank = %d\n\n',rank(Asym));
% --------------------------------------------------------
% High-precision eigenvalues
% --------------------------------------------------------
A = vpa(Asym,50);
evHP = eig(A);
% Convert to double for reliable sorting/classification
ev = double(evHP);
% Sort by real part
[~,idx] = sort(real(ev));
ev = ev(idx);
fprintf('Eigenvalues:\n\n');
for j = 1:length(ev)
fprintf('lambda_%d = %+.16e  %+.16e i\n', ...
j,real(ev(j)),imag(ev(j)));
end
% --------------------------------------------------------
% Count eigenvalues according to REAL PART
% --------------------------------------------------------
nplus  = nnz(real(ev) >  tol);
nminus = nnz(real(ev) < -tol);
nimag  = nnz(abs(real(ev)) <= tol);
fprintf('\n----------------------------------------\n');
fprintf('Classification by Re(lambda):\n');
fprintf('----------------------------------------\n');
fprintf('Right half-plane  Re(lambda) > 0 : %d\n',nplus);
fprintf('Left half-plane   Re(lambda) < 0 : %d\n',nminus);
fprintf('Imaginary axis    Re(lambda) = 0 : %d\n',nimag);
fprintf('TOTAL                           : %d\n', ...
nplus+nminus+nimag);
% --------------------------------------------------------
% Print only RHP eigenvalues
% --------------------------------------------------------
fprintf('\nEigenvalues in RIGHT half-plane:\n');
indRHP = find(real(ev) > tol);
if isempty(indRHP)
fprintf('None\n');
else
for j = indRHP.'
fprintf('%+.16e  %+.16e i\n', ...
real(ev(j)),imag(ev(j)));
end
end
% --------------------------------------------------------
% Print only LHP eigenvalues
% --------------------------------------------------------
fprintf('\nEigenvalues in LEFT half-plane:\n');
indLHP = find(real(ev) < -tol);
if isempty(indLHP)
fprintf('None\n');
else
for j = indLHP.'
fprintf('%+.16e  %+.16e i\n', ...
real(ev(j)),imag(ev(j)));
end
end
end
% ============================================================
% 3. Compare the two unstable counts
% ============================================================
A5     = double(vpa(subs(Gred,v6,5),50));
A40000 = double(vpa(subs(Gred,v6,40000),50));
ev5     = eig(A5);
ev40000 = eig(A40000);
nplus5 = nnz(real(ev5) > tol);
nplus40000 = nnz(real(ev40000) > tol);
fprintf('\n\n====================================================\n');
fprintf('FINAL COMPARISON\n');
fprintf('====================================================\n');
fprintf('Number of RHP eigenvalues at v6 = 5     : %d\n',nplus5);
fprintf('Number of RHP eigenvalues at v6 = 40000 : %d\n',nplus40000);
fprintf('\nExpected result:\n');
fprintf('n_+(Gred(5))     = 1\n');
fprintf('n_+(Gred(40000)) = 3\n');
fprintf('\n det Gred(v6) is nonzero for v6 in [5,40000],\n');

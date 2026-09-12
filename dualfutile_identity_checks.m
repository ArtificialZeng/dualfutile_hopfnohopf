function dualfutile_identity_checks()
%IDENTITY_CHECKS Exact symbolic checks of three identities used in the paper.

    clc;
    timer=tic;
    syms a b4 b10 c e1 e4 f7 f10 u2 u3 v5 v6 x8 x9 y11 y12 real

    %% 1. Full and reduced characteristic polynomials
    G=[ ...
        -a, 0, 0, -e1, 0, u2, 0, 0, y12;
        0, -b4-b10, 0, -e4, -f10, u3, v5, x9, y11;
        0, 0, -c, 0, -f7, 0, v6, x8, 0;
        -a, -b4, 0, -e1-e4, 0, u2+u3, v5+v6, 0, 0;
        0, -b10, -c, 0, -f7-f10, 0, 0, x8+x9, y11+y12;
        a, 0, 0, e1, 0, -u2-u3, 0, 0, 0;
        0, b4, 0, e4, 0, 0, -v5-v6, 0, 0;
        0, 0, c, 0, f7, 0, 0, -x8-x9, 0;
        0, b10, 0, 0, f10, 0, 0, 0, -y11-y12];

    Gred=[ ...
        0, 0, -u3, 0, y12, 0;
        0, 0, 0, v6, 0, -x9;
        a, 0, -a-e1-u2-u3, -e1, 0, 0;
        -b4, -b4, -e4, -b4-e4-v5-v6, -b4, 0;
        -b10, -b10, 0, -b10, -b10-f10-y11-y12, -f10;
        0, c, 0, 0, -f7, -c-f7-x8-x9];

    fullCoefficients=charpoly(G);
    reducedCoefficients=charpoly(Gred);
    assert(all(isAlways(fullCoefficients(1:7)==reducedCoefficients)) && ...
           all(isAlways(fullCoefficients(8:10)==0)), ...
        'The characteristic polynomials of G and Gred do not agree.');
    fprintf('1. det(lambda*I-G)=lambda^3 det(lambda*I-Gred): PASS\n');

    %% 2. Hurwitz identity
    syms c1 c2 c3 c4 c5 c6 real
    H=[ ...
        c1 c3 c5 0  0  0;
        1  c2 c4 c6 0  0;
        0  c1 c3 c5 0  0;
        0  1  c2 c4 c6 0;
        0  0  c1 c3 c5 0;
        0  0  1  c2 c4 c6];

    Delta2=expand(det(H(1:2,1:2)));
    Delta3=expand(det(H(1:3,1:3)));
    Delta4=expand(det(H(1:4,1:4)));
    Delta5=expand(det(H(1:5,1:5)));

    hurwitzIdentity=expand(Delta2*Delta5 ...
        -(c5*Delta2-c1^2*c6)*Delta4+c6*Delta3^2);
    assert(isAlways(hurwitzIdentity==0),'The Hurwitz identity failed.');
fprintf('2. Delta2*Delta5=(c5*Delta2-c1^2*c6)*Delta4-c6*Delta3^2: PASS\n');
    %% 3. P5,N5,P6,N6 rearrangement
    B1plus=expand(c3*Delta3-c1*c5*Delta2+c1^3*c6);
    assert(isAlways(expand(Delta5-c5*Delta4+c6*B1plus)==0), ...
        'The formula Delta5=c5*Delta4-c6*B1plus failed.');

    syms P5 N5 P6 N6 real
    decomposition=P6*Delta5-(N6*P5-N5*P6)*Delta4 ...
        -c6*(P5*Delta4-P6*B1plus);
    decomposition=expand(subs(decomposition,[c5 c6], ...
        [P5-N5 P6-N6]));
    assert(isAlways(decomposition==0),'The P5,N5,P6,N6 identity failed.');
fprintf('3. P6*Delta5=(N6*P5-N5*P6)*Delta4+c6*(P5*Delta4-P6*B1plus): PASS\n');
    fprintf('\nAll exact identity checks passed in %.1f seconds.\n',toc(timer));
end

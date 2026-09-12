function dualfutile_positivity_checks()
%PRELIMINARY_POSITIVITY_CHECKS
% Exact coefficientwise positivity checks used before the main positivity certificate.
%
% The script starts directly from the full 9-by-9 Jacobian G and constructs
% the first six coefficients of
%
%   det(lambda*I-G)
%      =lambda^3*(lambda^6+c1*lambda^5+...+c6)
%
% directly from the sixteen positive Jacobian parameters
%
%   a,b4,b10,c,e1,e4,f7,f10,u2,u3,v5,v6,x8,x9,y11,y12.
%
% It then verifies, by exact sparse-polynomial arithmetic, that:
%   (1) c0,c1,...,c4 are coefficientwise positive;
%   (2) Delta1,Delta2,Delta3 are coefficientwise positive;
%   (3) N6*P5-N5*P6 is coefficientwise positive;
%   (4) P5*Delta2-c1^2*P6 is coefficientwise positive;
%   (5) c1*c2-2*c3 is coefficientwise positive;
%   (6) the mass-action corrected numerator of Delta4 is
%       coefficientwise positive.
%
% The factor lambda^3 is taken from the conservation-law argument proved in
% the paper.  This script computes the first six coefficients directly from
% G, but does not separately construct c7,c8,c9 to verify that they vanish.
%
% No Symbolic Math Toolbox is used.  A polynomial is represented by packed
% uint64 exponent keys and integer-valued doubles.  Exactness is guarded
% before every potentially unsafe operation:
%   * before coefficient multiplication, the largest possible product is
%     asserted to be strictly smaller than flintmax;
%   * before accumarray summation, an absolute-sum bound is asserted to be
%     strictly smaller than flintmax, which bounds every partial sum even
%     when cancellation occurs;
%   * before packed-key addition, every resulting variable exponent is
%     asserted to be at most 15, preventing carry between 4-bit fields.

    clc;
    fprintf('Preliminary coefficientwise positivity checks\n\n');
    timer = tic;

    id.a=1;   id.b4=2;  id.b10=3; id.c=4;
    id.e1=5;  id.e4=6;  id.f7=7;  id.f10=8;
    id.u2=9;  id.u3=10; id.v5=11; id.v6=12;
    id.x8=13; id.x9=14; id.y11=15; id.y12=16;

    a=pvar(id.a);     b4=pvar(id.b4);   b10=pvar(id.b10); c=pvar(id.c);
    e1=pvar(id.e1);   e4=pvar(id.e4);   f7=pvar(id.f7);   f10=pvar(id.f10);
    u2=pvar(id.u2);   u3=pvar(id.u3);   v5=pvar(id.v5);   v6=pvar(id.v6);
    x8=pvar(id.x8);   x9=pvar(id.x9);   y11=pvar(id.y11); y12=pvar(id.y12);

    fprintf('Constructing the characteristic coefficients directly from G ...\n');

    % Species order: (A,B,C,E,F,U,V,X,Y).
    % We construct -G, since
    %   det(lambda*I-G)=det(lambda*I+(-G)).
    minusG=zero_matrix(9,9);

    minusG{1,1}=a;
    minusG{1,4}=e1;
    minusG{1,6}=pscale(u2,-1);
    minusG{1,9}=pscale(y12,-1);

    minusG{2,2}=padd(b4,b10);
    minusG{2,4}=e4;
    minusG{2,5}=f10;
    minusG{2,6}=pscale(u3,-1);
    minusG{2,7}=pscale(v5,-1);
    minusG{2,8}=pscale(x9,-1);
    minusG{2,9}=pscale(y11,-1);

    minusG{3,3}=c;
    minusG{3,5}=f7;
    minusG{3,7}=pscale(v6,-1);
    minusG{3,8}=pscale(x8,-1);

    minusG{4,1}=a;
    minusG{4,2}=b4;
    minusG{4,4}=padd(e1,e4);
    minusG{4,6}=pscale(padd(u2,u3),-1);
    minusG{4,7}=pscale(padd(v5,v6),-1);

    minusG{5,2}=b10;
    minusG{5,3}=c;
    minusG{5,5}=padd(f7,f10);
    minusG{5,8}=pscale(padd(x8,x9),-1);
    minusG{5,9}=pscale(padd(y11,y12),-1);

    minusG{6,1}=pscale(a,-1);
    minusG{6,4}=pscale(e1,-1);
    minusG{6,6}=padd(u2,u3);

    minusG{7,2}=pscale(b4,-1);
    minusG{7,4}=pscale(e4,-1);
    minusG{7,7}=padd(v5,v6);

    minusG{8,3}=pscale(c,-1);
    minusG{8,5}=pscale(f7,-1);
    minusG{8,8}=padd(x8,x9);

    minusG{9,2}=pscale(b10,-1);
    minusG{9,5}=pscale(f10,-1);
    minusG{9,9}=padd(y11,y12);

    ccoef=cell(1,7);
    ccoef{1}=pconst(1); % c0
    for k=1:6
        % The coefficient of lambda^(9-k) in det(lambda*I-G) is the
        % sum of the k-by-k principal minors of -G.  Since the full
        % polynomial contains the factor lambda^3, these are precisely
        % c1,...,c6 of its degree-six nonzero factor.
        ccoef{k+1}=principal_coefficient(minusG,k);
    end

    c0=ccoef{1}; c1=ccoef{2}; c2=ccoef{3}; c3=ccoef{4};
    c4=ccoef{5}; c5=ccoef{6}; c6=ccoef{7};

    expectedC=[1 16 97 280 394];
    for k=0:4
        assert_positive(sprintf('c_%d',k),ccoef{k+1},expectedC(k+1));
    end

    % Positive and negative parts of c5 and c6:
    %   c5=P5-N5, c6=P6-N6.
    P5=positive_part(c5); N5=negative_part(c5);
    P6=positive_part(c6); N6=negative_part(c6);
    assert_positive('P_5',P5,238);
    assert_positive('N_5',N5,2);
    assert_positive('P_6',P6,41);
    assert_positive('N_6',N6,9);
    fprintf('  c_5 split: %s positive and %s negative monomials.\n', ...
        comma_integer(numel(P5.k)),comma_integer(numel(N5.k)));
    fprintf('  c_6 split: %s positive and %s negative monomials.\n\n', ...
        comma_integer(numel(P6.k)),comma_integer(numel(N6.k)));

    fprintf('Constructing the first three Hurwitz determinants ...\n');
    Delta1=c1;
    Delta2=padd(pmul(c1,c2),c3,-1);
    Delta3=padd(pmul(c3,Delta2),pprod(c1,c1,c4),-1);
    Delta3=padd(Delta3,pmul(c1,c5));

    assert_positive('Delta_1',Delta1,16);
    assert_positive('Delta_2',Delta2,746);
    assert_positive('Delta_3',Delta3,41243);

    fprintf('\nConstructing the remaining preliminary quantities ...\n');
    ratioMargin=padd(pmul(N6,P5),pmul(N5,P6),-1);
    assert_positive('N_6 P_5-N_5 P_6',ratioMargin,1588);

    slopeMargin=padd(pmul(P5,Delta2),pprod(c1,c1,P6),-1);
    assert_positive('P_5 Delta_2-c_1^2 P_6',slopeMargin,72468);

    firstMargin=padd(pmul(c1,c2),c3,-2);
    assert_positive('c_1 c_2-2c_3',firstMargin,744);

    fprintf('\nConstructing Delta_4 and its mass-action correction ...\n');
    Delta30=padd(Delta3,pmul(c1,c5),-1);
    Delta4inner=padd(pmul(c2,Delta2),pmul(c1,c4),-2);

    Delta4=padd(pmul(c4,Delta30),pmul(c5,Delta4inner),-1);
    Delta4=padd(Delta4,pmul(c5,c5),-1);
    Delta4=padd(Delta4,pprod(c1,Delta2,c6));

    check_integer_poly(Delta4);
    assert(numel(Delta4.k)==1712913, ...
        'Unexpected support size for Delta_4.');
    assert(nnz(Delta4.c>0)==1712647 && nnz(Delta4.c<0)==266, ...
        'Unexpected sign pattern for Delta_4.');
    fprintf(['  Delta_4: %s positive and %s negative monomials ' ...
             '(expected raw expansion).\n'], ...
        comma_integer(nnz(Delta4.c>0)), ...
        comma_integer(nnz(Delta4.c<0)));

    PDelta4=positive_part(Delta4);
    negativeDelta4=negative_part(Delta4);

    % Every negative monomial of Delta4 contains a*c*f10*e4*v6*y12.
    negativeCore=pprod(a,c,f10,e4,v6,y12);
    coreKey=negativeCore.k;
    assert(numel(coreKey)==1 && all(divides_key(coreKey,negativeDelta4.k)), ...
    'A negative Delta_4 monomial does not contain a*c*f10*e4*v6*y12.');
    NDelta4=negativeDelta4;
    NDelta4.k=NDelta4.k-coreKey;
    check_integer_poly(NDelta4);

    % On the mass-action locus,
    %   f10*e4*v6*y12*(u2+u3)*(x8+x9)
    %      =e1*f7*u3*x9*(y11+y12)*(v5+v6).
    % Therefore the corrected numerator is
    %   (u2+u3)(x8+x9)P_Delta4
    %   -(y11+y12)(v5+v6)a*c*e1*f7*u3*x9*N_Delta4.
    left=pprod(padd(u2,u3),padd(x8,x9),PDelta4);
    right=pprod(padd(y11,y12),padd(v5,v6), ...
        a,c,e1,f7,u3,x9,NDelta4);
    correctedDelta4=padd(left,right,-1);

    assert_positive( ...
        '(u2+u3)(x8+x9)P_Delta4-(y11+y12)(v5+v6)ace1f7u3x9 N_Delta4', ...
        correctedDelta4,4295079);

    fprintf('\nAll requested coefficientwise positivity checks passed.\n');
    fprintf('Elapsed time: %.1f seconds.\n',toc(timer));
end

% ======================================================================
% Characteristic coefficients
% ======================================================================

function C = principal_coefficient(M,k)
    n=size(M,1);
    subsets=nchoosek(1:n,k);
    C=pzero();
    for row=1:size(subsets,1)
        ix=subsets(row,:);
        C=padd(C,pdet(M(ix,ix)));
    end
end

function D = pdet(M)
    n=size(M,1);
    permutations=perms(1:n);
    D=pzero();
    for row=1:size(permutations,1)
        permutation=permutations(row,:);
        term=pconst(permutation_sign(permutation));
        for k=1:n
            term=pmul(term,M{k,permutation(k)});
            if isempty(term.k), break; end
        end
        D=padd(D,term);
    end
end

function s = permutation_sign(p)
    inversions=0;
    for i=1:numel(p)
        for j=i+1:numel(p)
            inversions=inversions+(p(i)>p(j));
        end
    end
    s=1-2*mod(inversions,2);
end

% ======================================================================
% Positivity and support checks
% ======================================================================

function assert_positive(name,P,expectedSupport)
    check_integer_poly(P);
    assert(~isempty(P.k),sprintf('%s is the zero polynomial.',name));
    assert(all(P.c>0),sprintf('%s is not coefficientwise positive.',name));
    if nargin>=3
        assert(numel(P.k)==expectedSupport, ...
            sprintf('Unexpected support size for %s.',name));
    end
    fprintf('  %-58s PASS  (%s positive monomials)\n', ...
        name,comma_integer(numel(P.k)));
end

function P = positive_part(P)
    keep=P.c>0;
    P.k=P.k(keep);
    P.c=P.c(keep);
end

function P = negative_part(P)
    keep=P.c<0;
    P.k=P.k(keep);
    P.c=-P.c(keep);
end

function yes = divides_key(divisor,keys)
    yes=true(size(keys));
    for variable=1:16
        yes=yes & getexp(keys,variable)>=getexp(divisor,variable);
    end
end

% ======================================================================
% Exact sparse-polynomial arithmetic
% ======================================================================

function P = pzero()
    P=struct('k',zeros(0,1,'uint64'),'c',zeros(0,1));
end

function P = pconst(c)
    if c==0
        P=pzero();
    else
        P=struct('k',uint64(0),'c',double(c));
    end
end

function P = pvar(i)
    P=struct('k',bitshift(uint64(1),4*(i-1)),'c',1);
end

function P = pscale(P,s)
    if s==0
        P=pzero();
        return
    end
    assert_exact_product_bound(P.c,s);
    P.c=P.c*s;
end

function C = padd(A,B,scaleB)
    if nargin<3, scaleB=1; end
    if isempty(A.k), C=pscale(B,scaleB); return; end
    if isempty(B.k), C=A; return; end
    scaledB=pscale(B,scaleB);
    C=pnormalize([A.k;scaledB.k],[A.c;scaledB.c]);
end

function C = psum(varargin)
    C=pzero();
    for k=1:nargin
        C=padd(C,varargin{k});
    end
end

function C = pprod(varargin)
    C=pconst(1);
    sizes=cellfun(@(q)numel(q.k),varargin);
    [~,order]=sort(sizes);
    for k=order
        C=pmul(C,varargin{k});
    end
end

function C = pmul(A,B)
    if isempty(A.k) || isempty(B.k)
        C=pzero();
        return
    end
    if numel(A.k)>numel(B.k)
        T=A; A=B; B=T;
    end

    % A packed monomial uses four bits for each of sixteen exponents.
    % Checking each field before adding keys prevents both field carry and
    % uint64 overflow.  Thus every exponent produced here lies in [0,15].
    assert_product_exponent_bound(A.k,B.k);
    assert_exact_product_bound(A.c,B.c);

    if numel(A.k)==1
        C.k=B.k+A.k;
        C.c=B.c*A.c;
        check_integer_poly(C);
        return
    end

    pairBudget=1.5e6;
    numberB=numel(B.k);
    block=max(1,floor(pairBudget/numberB));
    numberChunks=ceil(numel(A.k)/block);
    chunks=cell(numberChunks,1);
    chunk=0;

    for first=1:block:numel(A.k)
        chunk=chunk+1;
        ii=first:min(first+block-1,numel(A.k));
        keys=A.k(ii)+B.k.';
        coefficients=A.c(ii)*B.c.';
        chunks{chunk}=pnormalize(keys(:),coefficients(:));
    end

    while numel(chunks)>1
        next=cell(ceil(numel(chunks)/2),1);
        out=0;
        for k=1:2:numel(chunks)
            out=out+1;
            if k==numel(chunks)
                next{out}=chunks{k};
            else
                next{out}=padd(chunks{k},chunks{k+1});
            end
        end
        chunks=next;
    end
    C=chunks{1};
end

function P = pnormalize(keys,coefficients)
    keep=coefficients~=0;
    keys=keys(keep);
    coefficients=coefficients(keep);
    if isempty(keys), P=pzero(); return; end

    % If the sum of the absolute values of all input coefficients is below
    % flintmax, then every partial sum used for any monomial group is below
    % flintmax as well, independently of summation order and cancellation.
    assert_exact_sum_bound(coefficients);

    [keys,order]=sort(keys);
    coefficients=coefficients(order);
    first=[true;keys(2:end)~=keys(1:end-1)];
    assert(numel(first)<flintmax,'The grouping index exceeds flintmax.');
    group=cumsum(first);
    coefficients=accumarray(group,coefficients,[],@sum);
    keys=keys(first);
    keep=coefficients~=0;
    P=struct('k',keys(keep),'c',coefficients(keep));
    check_integer_poly(P);
end

function check_integer_poly(P)
    assert(all(P.c==fix(P.c)) && all(abs(P.c)<flintmax), ...
        'A coefficient is no longer an exactly represented integer.');
end

function assert_exact_product_bound(A,B)
    % All inputs have already been created by guarded exact operations.
    assert(all(A(:)==fix(A(:))) && all(B(:)==fix(B(:))) && ...
           all(abs(A(:))<flintmax) && all(abs(B(:))<flintmax), ...
        'A multiplication input is outside the exact-integer range.');
    if isempty(A) || isempty(B), return; end

    maxA=max(abs(A(:)));
    maxB=max(abs(B(:)));
    limit=uint64(flintmax)-uint64(1);
    intA=uint64(maxA);
    intB=uint64(maxB);
    assert(double(intA)==maxA && double(intB)==maxB, ...
        'A multiplication input is outside the exact-integer range.');
    if intA~=0
        assert(intB<=idivide(limit,intA,'floor'), ...
            'An intermediate coefficient product may reach flintmax.');
    end
end

function assert_exact_sum_bound(coefficients)
    assert(all(coefficients(:)==fix(coefficients(:))) && ...
           all(abs(coefficients(:))<flintmax), ...
        'A summation input is outside the exact-integer range.');
    if isempty(coefficients), return; end

    maxCoefficient=max(abs(coefficients(:)));
    integerMaximum=uint64(maxCoefficient);
    assert(double(integerMaximum)==maxCoefficient, ...
        'A summation input is outside the exact-integer range.');
    if integerMaximum==0, return; end

    % numel(coefficients)*max(abs(coefficients)) is an upper bound for the
    % absolute sum of every possible subset, hence for every intermediate
    % partial sum formed by accumarray.  Integer division avoids evaluating
    % that product in floating-point arithmetic.
    limit=uint64(flintmax)-uint64(1);
    allowedCount=idivide(limit,integerMaximum,'floor');
    assert(uint64(numel(coefficients))<=allowedCount, ...
        'An intermediate coefficient sum may reach flintmax.');
end

function assert_product_exponent_bound(keysA,keysB)
    for variable=1:16
        maximumA=max(double(getexp(keysA,variable)));
        maximumB=max(double(getexp(keysB,variable)));
        assert(maximumA+maximumB<=15, ...
            sprintf('The packed exponent field for variable %d may overflow.', ...
                    variable));
    end
end

function M = zero_matrix(m,n)
    M=cell(m,n);
    for i=1:m
        for j=1:n
            M{i,j}=pzero();
        end
    end
end

function e = getexp(keys,i)
    e=uint8(bitand(bitshift(keys,-4*(i-1)),uint64(15)));
end

function text = comma_integer(n)
    text=sprintf('%.0f',n);
    first=mod(numel(text),3);
    if first==0, first=3; end
    parts={text(1:first)};
    for k=first+1:3:numel(text)
        parts{end+1}=text(k:k+2); %#ok<AGROW>
    end
    text=strjoin(parts,',');
end

function [S,G, RSE] = SSNMTF(R, A, k, gamma, max_iter, random_init, G0,FOUT, verbose)
% Function for solving Graph Regularized Symmetric Non-Negative Matrix Tri-Factorization
% -------------------------------------------------------------------------------------------------------------
%
% Noel Malod-Dognin, University College London.
%
% based on the function factorization written by
% Janez Povh and Vladimir Gligorijevic
% --------------------------------------------------------------------------------------------------------------
%
% \min[ \sum_i( ||R_i - G*S_i*G^T||^2 ) + \sum_j( \gamma_j*Trace(G^T*LA_j*G) ] such that G\ge 0
%
% [Input]:
%     R: <1D cell array> of i n x n symmetric relational matrix on a set of n objects,
%     A: <1D cell array> of j n x n adjacency matrices connecting the n objects
%     k: <int>, rank parameters for matrices S_i (S_i are k x k symmetric matrices
%     gamma: <1 x j float> values for parameter gamma (one per adjacency matrix in A)
%     max_iter: <int>, predefined number of iterations
%     G0 ... n x k matrix - initial matrix
%     FOUT if of output to save the results
%
% [Output]:
%     S: <1D Cell>, k(node types) x k(node types) blocks, compressed matrix (e.g., S{i} = Si, k x k)
%     G: r x r nonnegative matrix shared over all R_i
%     RSE: double   RSE = \sum_i||R_i - G*S_i*G^T||^2/\sum_i||R_i||^2
% --------------------------------------------------------------------------------------------------------------

if nargin <= 6
    FOUT = 0;
end;

if nargin <= 7
    verbose = 0;
end;

% order of matrices R_i
nR=[];
r=length(R);
for i=1:length(R)
    nR=[nR;size(R{i})];
end;
if min(nR(:))~=max(nR(:))
    error('Matrices Ri should be of the same order');
else
    n=min(nR(:));
end;

% computing norm of R matrix
norm_R = 0;
Rsum=zeros(n);
for i=1:r
  norm_R = norm_R + (norm(R{i},'fro'))^2;
  Rsum = Rsum + R{i};
end;

% Random initialization of G matrix
if random_init
% fprintf('-Initializations of G matrix....\n');
    G = rand(n,k);
end;
else
% fprintf('-Initialization of G matrix finished!\n\n');
%

% SVD based initialization of G, from
% [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
%     start for nonnegative matrix factorization, Pattern Recognition,
%     Elsevier

    if nargin==6
        G=G0;
    else
        [U,D]=eig(Rsum);
        diagDpos=diag(abs(D));
        [Dsort,ind]=sort(diagDpos,1,'descend');%

        trD=sum(diagDpos);
        trDi=0;
        K=1;
        for i=1:length(Dsort)
            trDi = trDi + Dsort(i);
            if trDi/trD > 0.9
                K=i
                break;
            end;
        end;
        k=min(k,K);
        Ginit=[];
        for i=1:k
            xx=U(:,ind(i))*diagDpos(ind(i));
            xp=max(xx,0);
            xn=xp-xx;
            if norm(xp)>norm(xn)
              Ginit = [Ginit xp];
            else
              Ginit = [Ginit xn];
            end
        end;
        Gzeros=Ginit<1e-10;
        G=Ginit+1e-10*Gzeros;
    end;
end;

%Graph Regularization
a=length(A);
for ii=1:a
    L_pos{ii} = diag(sum(A{ii},1));
    L_neg{ii} = A{ii};
end;


J_old = 0; %initialization
%Iterations
RSE = 1;

if verbose
    fprintf('| Iteration | Delta_R |  RSE | Rel_Var_Cost | KKT_err | Time \n');
end;
if FOUT
    fprintf(FOUT,'| Iteration | Delta_R | RSE | Rel_Var_Cost | KKT_err | Time \n');
end;
tic;


for iter=1:max_iter

    GtG = G'*G + eps;
    GtG_inv = inv(GtG);

    % Update S
    for i=1:r
        S{i} = GtG_inv*G'*R{i}*G*GtG_inv;
    end;

    % Initialize sum (numerator and denominator)
    Ge = zeros(n,k);
    Gd = zeros(n,k);

    % Update G
    for i=1:r
      RiGSi = R{i}*G*S{i};
      RiGSi_pos = (abs(RiGSi)+RiGSi)/2.0;
      RiGSi_neg = (abs(RiGSi)-RiGSi)/2.0;

      SiGtGSi = S{i}*GtG*S{i};
      SiGtGSi_pos = (abs(SiGtGSi)+SiGtGSi)/2;
      SiGtGSi_neg = (abs(SiGtGSi)-SiGtGSi)/2;

      GSiGtGSi_pos = G*SiGtGSi_pos;
      GSiGtGSi_neg = G*SiGtGSi_neg;

      Ge = Ge + RiGSi_pos + GSiGtGSi_neg;
      Gd = Gd + RiGSi_neg + GSiGtGSi_pos+eps*ones(n,k);
    end;

     % Adding constraints and computing new values for Gi
    for t=1:a
        Ge = Ge + gamma(t)*L_neg{t}*G;
        Gd = Gd + gamma(t)*L_pos{t}*G;
    end;

    % Updating G
    G=G.*sqrt(Ge./Gd);

    % Computing the relative square error (RSE) every 10th iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (n-1)th iteration
    if mod(iter+1,10) == 0

        relation_error = 0;
        penalty = 0;
        for i=1:a
            penalty = penalty + gamma(a)*trace(G'*(L_pos{i}-L_neg{i})*G);
        end;
        for i=1:r
            relation_error = relation_error + (norm( R{i} - G*S{i}*G', 'fro' ))^2;
        end;
        %fprintf('REL = %f, PEN = %f: ', relation_error, penalty)
        % objective (cost) function
        J_old = relation_error + penalty;

    end;

    % nth iteration
    if mod(iter,10) == 0

        relation_error = 0;
        penalty = 0;
        for i=1:a
            penalty = penalty + gamma(a)*trace(G'*(L_pos{i}-L_neg{i})*G);
        end;
        for i=1:r
           relation_error = relation_error + (norm( R{i} - G*S{i}*G', 'fro' ))^2;
        end;
        % objective (cost) function
        J_new = relation_error + penalty;


        % Errors
        RSE = relation_error/norm_R;
        rel_var = abs(J_new - J_old)/abs(J_old);

        % KKT error - how KKT conditiona are violated
        Err1=zeros(size(G));
        Err2=zeros(size(S{1}));
        GtG=G'*G;
        for iix=1:r
           Err1=Err1+R{iix}*G*S{iix}-G*S{iix}*GtG*S{iix};
           Err2=Err2+G'*R{iix}*G-GtG*S{iix}*GtG;
        end;
        Err3=Err1(:)'*G(:);
        KKT_err=max(max(norm(Err1,'fro'),norm(Err2,'fro')),Err3);
        % Writing output
        ttoc=toc;
        if verbose
            fprintf('%d %0.5e %0.5e %0.5e %0.5e %0.5e\n', iter,  J_new, RSE, rel_var,KKT_err,ttoc);
        end;
        if FOUT
          fprintf(FOUT,'%d %0.5e %0.5e %0.5e %0.5e %0.5e\n', iter,  J_new, RSE, rel_var,KKT_err,ttoc);
        end;


        if (RSE <= 1e-2 | rel_var <= 1e-8) % checking for convergence
            break;
        end;

    end;

end;
%keyboard;

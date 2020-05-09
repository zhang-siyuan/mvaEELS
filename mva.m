%MVA 
% mva is a Matlab script that performs PCA and/or NNMF analysis with 
% pre-defined or user-defined numbers of components through an interactive
% interface. 
% Author: Siyuan Zhang (<a
% href="mailto:siyuan.zhang@mpie.de">siyuan.zhang@mpie.de</a>)
% Please cite this paper: https://doi.org/10.1093/jmicro/dfx091

function [spect,weight,x,y,scree,residue,A] = mva(A,e,x,y,norm,nComp)
if nargin < 6
    nComp = 0;
end
residue = [];
scree = [];
set(0,'DefaultFigureColormap',hot);
A = double(A');
A = A - min(0,min(A(:)));
A0 = A;
Nx = length(x);
Ny = length(y);
Ne = length(e);
if nComp <= 0 % nComp > 0: no PCA; = 0: PCA&choose; = -1: PCA&NMF; = -2: PCA only;
    sOption = -nComp;
%% PCA
screenSize = get(0,'ScreenSize');
figID = figure('KeyPressFcn',@(obj,evt) 0,'Position',[1 screenSize(4)/2.5 screenSize(3) screenSize(4)/2]);
key = 13;
maxComp = 1;
while maxComp
    figure(figID)
    if key == 13
        if (nargin >= 5) && norm
            A1 = sqrt(abs(mean(A,1)));
            A2 = sqrt(abs(mean(A,2)));
            a1 = 1./A1;
            a1(a1 == inf) = 0;
            a2 = 1./A2;
            a2(a2 == inf) = 0;
        else
            A1 = ones(1,size(A,2));
            A2 = ones(size(A,1),1);
            a1 = A1;
            a2 = A2;
        end
        A = sparse(diag(a2))*A*sparse(diag(a1));
        A(isnan(A)) = 0;
        pix = Nx*Ny;
        [wA,scree,sA] = svd(A);
        scree = diag(scree);
        maxComp = length(scree);
        for iComp = 1:maxComp
            if sum(wA(:,iComp)) < 0
                wA(:,iComp) = -wA(:,iComp);
                sA(:,iComp) = -sA(:,iComp);
            end
        end
        A = sparse(diag(A2))*A*sparse(diag(A1));
        wA = sparse(diag(A2))*wA;
        sA = sparse(diag(A1))*sA;
        nComp = 1;
        xRange = 1:Nx;
        yRange = 1:Ny;
    end
    subplot(1,3,1)
    loglog(nComp:maxComp,scree(nComp:maxComp),'o')
    axis tight
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    title('PCA Scree Plot');    
%     title([sprintf('PCA Scree Plot\nInclude first %d components\n',nComp),...
%         '[\leftarrow],[\rightarrow]:Navigate to other components',sprintf('\n'),...
%         '[\uparrow],[\downarrow]:Adjust last component in view, [Q]:Proceed']);
    subplot(1,3,2)
    wa = reshape(wA,Nx,Ny,pix);
    imagesc(x(xRange),y(yRange),wa(xRange,yRange,nComp)');
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    axis equal tight off
    title(sprintf('Weighting of component %d\n',nComp));
%     title([sprintf('Weighting of component %d\n',nComp),...
%         sprintf('[W][A][S][D]:Cut dimensions%d * %d\n',length(xRange),length(yRange)),...
%         '[R]:Reset cuts, [Enter]:Accept cuts',sprintf('\n'),...
%         '[B][X][Y]:Bin dimensions']);
    subplot(1,3,3)
    plot(e,sA(:,nComp));
    axis tight
    set(gca,'FontSize',12)
    set(gca,'lineWidth',2)
    title(sprintf('Spectrum'));
    waitfor(gcf,'CurrentCharacter');
    key = uint8(get(gcf,'CurrentCharacter'));
    set(gcf,'CurrentCharacter',char(0))
    if key == 13 % Enter: compute PCA
        A = reshape(A,Nx,Ny,Ne);
        A = A(xRange,yRange,:);
        x = x(xRange);
        y = y(yRange);
        Nx = length(xRange);
        Ny = length(yRange);
        A = reshape(A,Nx*Ny,Ne);
    elseif key == 28 % Left arrow
        nComp = max(1,nComp-1);
    elseif key == 29 % Right arrow
        nComp = min(maxComp,nComp+1);
    elseif key == 30 % Up arrow
        maxComp = min(maxComp*2,length(scree));
    elseif key == 31 % Down arrow
        maxComp = ceil(max(maxComp/2,nComp));
    elseif key == 97 % A: Left crop
        if length(xRange) > 1
            xRange(1) = [];
        end
    elseif key == 100 % D: Right crop
        if length(xRange) > 1
            xRange(end) = [];
        end
    elseif key == 119 % W: Up crop
        if length(yRange) > 1
            yRange(1) = [];
        end
    elseif key == 115 % S: Down crop
        if length(yRange) > 1
            yRange(end) = [];
        end
    elseif key == 114 % R: reset
        xRange = 1:Nx;
        yRange = 1:Ny;
    elseif key == 113 % Q: quit
        maxComp = 0;
    elseif key == 98 % B: bin by 2
        [A,Nx,Ny] = binCube(A,Nx,Ny,2);
        x = x(2:2:2*Nx);
        y = y(2:2:2*Ny);
        key = 13;
    elseif key == 120 % X: bin by 2
        [A,Nx,Ny] = binCube(A,Nx,Ny,[2,1]);
        x = x(2:2:2*Nx);
        key = 13;
    elseif key == 121 % Y: bin by 2
        [A,Nx,Ny] = binCube(A,Nx,Ny,[1,2]);
        y = y(2:2:2*Ny);
        key = 13;
    elseif key == 127 % Del: component removal
%         A = wA(:,[1:nComp-1,nComp+1:end])*diag(scree([1:nComp-1,nComp+1:end]))*sA([1:nComp-1,nComp+1:end],:);
        A = A - wA(:,nComp)*scree(nComp)*sA(:,nComp)';
        key = 13;
    end
end
close(figID)
else
    sOption = 1;
end
if sOption == 0
    sOption = menu('Continue with this option',{'NNMF analysis','Quit'});
end
if sOption == 2
    spect = sA(:,1:nComp);
    weight = diag(scree(1:nComp))*wA(:,1:nComp)';
    A = spect*weight;
    A = floor(A);
    return
end
%% NNMF
if (nargin >= 5) && norm
    A1 = sqrt(mean(A,1));
    A2 = sqrt(mean(A,2));
    a1 = 1./A1;
    a1(a1 == inf) = 0;
    a2 = 1./A2;
    a2(a2 == inf) = 0;
else
    A1 = ones(1,size(A,2));
    A2 = ones(size(A,1),1);
    a1 = A1;
    a2 = A2;
end
A = sparse(diag(a2))*A*sparse(diag(a1));
A(isnan(A)) = 0;
A(A<0) = 0;
[weight,spect,residue.R] = nnmf_serial(A,nComp,20,20,1e-10);
A = sparse(diag(A2))*A*sparse(diag(A1));
weight = sparse(diag(A2))*weight;
spect = spect*sparse(diag(A1)); %%double check
for iComp = 1:nComp
    scaleF = sum(weight(:,iComp));%/Nx/Ny;
    spect(iComp,:) = spect(iComp,:)*scaleF;
    weight(:,iComp) = weight(:,iComp)/scaleF;
end
%% Component selection
useComp = ones(1,nComp+1);
screenSize = get(0,'ScreenSize');
figID = figure('KeyPressFcn',@(obj,evt) 0,'Position',[1 screenSize(4)/2.5 screenSize(3) screenSize(4)/2]);
iSelect = 1;
while iSelect
    figure(figID)
    if iSelect < nComp+1
        subplot(1,2,1)
        plot(e,spect(iSelect,:),'lineWidth',2)
        axis tight
        set(gca,'FontSize',12)
        set(gca,'lineWidth',2)
        if useComp(iSelect)
            tempStr1 = '';
            tempStr2 = 'Disregard';
        else
            tempStr1 = 'not ';
            tempStr2 = 'Select';
        end
        title([sprintf('Spectral Component %d: ',iSelect),tempStr1,...
            sprintf('selected\n'),'[Enter]:',tempStr2,' the component',...
            sprintf('\n'),'[\leftarrow],[\rightarrow]:Navigate to other components',...
            sprintf('\n'),'[Q]:Proceed']);
        subplot(1,2,2)
        imagesc(x,y,reshape(weight(:,iSelect),Nx,Ny)')
        axis equal tight
        set(gca,'FontSize',12)
        set(gca,'lineWidth',2)
        title('Weighting map');
    else            
        spectComp = find(useComp(1:nComp) > 0);
        subplot(1,2,1)
        imagesc(e,1:Nx*Ny,weight(:,spectComp)*spect(spectComp,:));
        title(['Sum components: ',num2str(spectComp),sprintf('\n'),...
            '[\leftarrow]:Navigate to other components, [Q]:Proceed']);
        colorbar
        subplot(1,2,2)
        residue.data = weight(:,spectComp)*spect(spectComp,:)-A;
        imagesc(e,1:Nx*Ny,residue.data);
        title('Residual from original spectrum image');
        colorbar
    end
    waitfor(gcf,'CurrentCharacter');
    key = uint8(get(gcf,'CurrentCharacter'));
    set(gcf,'CurrentCharacter',char(0))
    if key == 13 % Enter: select/deselect
        useComp(iSelect) = ~useComp(iSelect);
    elseif key == 28 % Left arrow
        iSelect = max(iSelect-1,1);
    elseif key == 29 % Right arrow
        iSelect = min(iSelect+1,nComp+1);
    elseif key == 30 % Up arrow
        iSelect = 1;
    elseif key == 31 % Down arrow
        iSelect = nComp;
    elseif key == 113 % Q: quit
        if sum(useComp) > 0
            iSelect = 0;
        else
            disp('Warning: select at least 1 component');
        end
    elseif key == 115 % S: sum
        iSelect = nComp+1;
    end
end
close(figID)
useComp(nComp+1) = [];
spectComp = find(useComp(1:nComp) > 0);
residueComp = find(useComp(1:nComp) == 0);
for iComp = residueComp
    residue.(['Spect',num2str(iComp)]) = spect(iComp,:);
    residue.(['Weight',num2str(iComp)]) = reshape(weight(:,iComp),Nx,Ny)';
end
spect = spect(spectComp,:)';
weight = weight(:,spectComp)';
A = floor(A');
end

function [output,nx,ny] = binCube(data,Nx,Ny,bxy)
if length(bxy) == 1
    bx = bxy;
    by = bxy;
else
    bx = bxy(1);
    by = bxy(2);
end
Nz = size(data,2);
nx = Nx/bx;
rx = 0;
if nx ~= floor(nx)
    nx = floor(nx);
    rx = Nx - nx*bx;
end
ny = Ny/by;
ry = 0;
if ny ~= floor(ny)
    ny = floor(ny);
    ry = Ny - ny*by;
end
data = reshape(data,Nx,Ny,Nz);
data(end+1-rx:end,:,:) = [];
data(:,end+1-ry:end,:) = [];
data = reshape(data,bx,nx,by,ny,Nz);
output = reshape(sum(sum(data,1),3),[nx*ny Nz]);
end

function [W,H,D] = nnmf_serial(a,k,init_iter,max_iter,tol)
%NORM_NNMF_SERIAL Normalised non-negative matrix factorisation (serial)
% [W,H,D] = NORM_NNMF_SERIAL(A,K,MAX_ITER,TOL) normalises the input matrix 
% A, and decomposes into K components by NNMF. The best of INIT_ITER 
% random seeds is fed into a SK_LEARN NNMF algorithm until the maximum 
% number of iterations MAX_ITER or the tolerance TOL is reached.
% Authors: Wai Yuen Fu (<a
% href="mailto:wyf22@cam.ac.uk">wyf22@cam.ac.uk</a>), Siyuan Zhang (<a
% href="mailto:siyuan.zhang@mpie.de">siyuan.zhang@mpie.de</a>)
% Copyright 2014-2015 University of Cambridge
%
% See also NNMF.

if (nargin <= 2)
    init_iter = 10;
end
if (nargin <= 3)
    max_iter = 100;
end
if (nargin <= 4)
    tol = 1e-6;
end

a = double(a);

minv = min(min(a));
if (minv < 0)
    a = a-minv;
end

% aGh = sqrt(mean(a,2));
% aGi = 1./aGh;
% aGi(aGi == inf) = 0;
% bHh = sqrt(mean(a,1));
% bHi = 1./bHh;
% bHi(bHi == inf) = 0;
% 
% aw = sparse(diag(aGi))*a*sparse(diag(bHi));
% aw(isnan(aw)) = 0;

[W,H,D] = nnmf_sk(a,k,init_iter,max_iter,tol);

% W = sparse(diag(aGh))*W;
% H = H*sparse(diag(bHh));

fprintf('Final root mean square residual = %.6f\n',D);

end

function [W,H,D] = init_nnmf(a,k,tol,max_iter)   
    D = 1e10;
    for q = 1:max_iter
        [W0,H0] = init_nnmf_one(a,k);
        [W1,H1,D1] = nnmf(a,k,'w0',W0,'h0',H0,'algorithm','als','options',...
            statset('TolFun',tol*1000,'Display','off'));
        if (D1 < D)
            W = W1;
            H = H1;
            D = D1;
        end
    end
    fprintf(' Seed root mean square residual = %.6f\n',D);
    
end

function [W,H] = init_nnmf_one(a,k,tol,variant)
if nargin <= 2
    tol = 1e-6;
end
if nargin <= 3
    variant = 0;
end
    [U,S,V] = rsvd(a,k);
    V = V';
    S = diag(S);
    W = zeros(size(U));
    H = zeros(size(V));
    W(:,1) = sqrt(S(1))*abs(U(:,1));
    H(1,:) = sqrt(S(1))*abs(V(1,:));
    
    for j = 2:k
        x = U(:,j); x_p = x; x_n = -x;
        y = V(j,:); y_p = y; y_n = -y;
        % extract positive and negative parts of column vectors
        x_p(x_p<0) = 0; x_n(x_n<0) = 0;
        y_p(y_p<0) = 0; y_n(y_n<0) = 0;
        % and their norms
        x_p_nrm = norm(x_p); x_n_nrm = norm(x_n);
        y_p_nrm = norm(y_p); y_n_nrm = norm(y_n);
        m_p = x_p_nrm * y_p_nrm;
        m_n = x_n_nrm * y_n_nrm;
        % choose update
        if (m_p > m_n) 
            u = x_p / x_p_nrm;
            v = y_p / y_p_nrm;
            sigma = m_p;
        else
            u = x_n / x_n_nrm;
            v = y_n / y_n_nrm;
            sigma = m_n;
        end
        lbd = sqrt(S(j) * sigma);
        W(:, j) = lbd * u;
        H(j, :) = lbd * v;
        
    end
    
    W(W<tol) = 0;
    H(H<tol) = 0;
    
    if variant == 1
        avg = mean(mean(a));
        W(W == 0) = avg;
        H(H == 0) = avg;
    elseif variant == 2
        avg = mean(mean(a));
        W(W == 0) = abs(avg * randn(len(W(W == 0))) / 100);
        H(H == 0) = abs(avg * randn(len(H(H == 0))) / 100);
    end
end

function [W,H,D] = nnmf_sk(a,k,init_iter,max_iter,tol)
tol_inner = tol/10;
% SKLearn Version
[W0,H0,D0] = init_nnmf(a,k,tol,init_iter);

for q = 1:max_iter
    [W,H,D] = nnmf(a,k,'w0',W0,'h0',H0,'algorithm','mult','options',...
        statset('MaxIter',max_iter,'Display','off','TolFun',tol_inner,'TolX',tol_inner));
    if (D < D0)
        W0 = W;
        H0 = H;
        
        if (D+tol >= D0)
            break;
        else
            D0 = D;
        end
    else
        W = W0;
        H = H0;
        D = D0;
        break;
    end
    
end
end

function [U,S,V] = rsvd(A,K)
%-------------------------------------------------------------------------------------
% random SVD
% Extremely fast computation of the truncated Singular Value Decomposition, using
% randomized algorithms as described in Halko et al. 'finding structure with randomness
%
% usage : 
%
%  input:
%  * A : matrix whose SVD we want
%  * K : number of components to keep
%
%  output:
%  * U,S,V : classical output as the builtin svd matlab function
%-------------------------------------------------------------------------------------
% Antoine Liutkus  (c) Inria 2014

[~,N] = size(A);
P = min(2*K,N);
X = randn(N,P);
Y = A*X;
W1 = orth(Y);
B = W1'*A;
[W2,S,V] = svd(B,'econ');
U = W1*W2;
K=min(K,size(U,2));
U = U(:,1:K);
S = S(1:K,1:K);
V=V(:,1:K);
end
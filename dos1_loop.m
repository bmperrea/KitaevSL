function [En,Ev,dd,Dd,Ixx,Ixy,Ixy2] = dos1_loop(T,H,Rxx,Rxy,RxyA,RxyB,bins,emax)
%Computes the DOS and three Raman spectra for the given Hamiltonian,
%temperatures array T, and Raman operators Rs, all in the form of submatrices for
%a case assuming sublattice symmetry (so H is really D from my notes)

%Inputs are temperature T in units of J (kitaev exchange) and a sparse
%matrix for the Hamiltonian H and a cell of sparse matrices for the Raman
%operators Rs with {Rxx,Rxy,R[xy]A,R[xy]B} and the last two are for the two
%sublattices

%dd   DOS
%Dd   two-particle DOS

%Hard coded data
%bins = 200;
%emax = 12; 

%useful data
S = size(H,1);
NN = 2*S;
Ev = (1:bins)'*emax/bins;
Z = NN*(emax/2)/bins;

%toc
%tic

%memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We play some tricks to keep the memory down to a minimum.

%Diagonalizing
H2 = full(H);
clearvars H

% H = U*D*V'   ==>    D = U'*H*V
[U,D,V] = svd(H2);
%memory
clearvars H2

%toc
%tic

% en2 = flip(diag(D));
% clearvars D
% U = fliplr(U);
% V = fliplr(V);

en2 = 2*diag(D);
clearvars D

En = -sum(en2)/2;

%Store the energies in a big matrix
ent0 = bsxfun(@plus,en2',en2);%repmat(en2,1,S) + repmat(en2',S,1);
%ins0 = (emax/2>=ent0)&(ent0>=0);
%ent1 = bsxfun(@minus,en2',en2);%-repmat(en2,1,S) + repmat(en2',S,1);
%ins1 = (emax/2>=ent1)&(ent1>=0);

%Store the 1p-DOS as well
%ins2 = (emax/4>=en2) & en2>=0;
histv0 = histwv_loop(en2(:),ones(size(en2(:))),0,emax/2,bins);
dd = 2*histv0/Z;

W0 = zeros(S^2,numel(T),'single');

if numel(T) == 1 && T == 0
    
    Rxx2 = U'*Rxx*V;
    R0 = -Rxx2 + Rxx2';
 %   R1 = -Rxx2 - Rxx2';
    clearvars Rxx2 Rxx
    % See Eq. 21 of 1507.01639
    R0 = pi*real( conj(R0).*R0 );
%    R1 = 2*pi*real( R1.*(R1.') );
%    R02 = bsxfun(@times,bsxfun(@times,R0,nf),nf');
    [histw0, histv0, subs0, inds0] = ...
        histwv_loop(reshape(ent0,S^2,1),reshape(R0,S^2,1),0,emax,bins);
%    R12 = bsxfun(@times,bsxfun(@times,R1,nf),1-nf');
 %   [histw2, histv2] = histwv_loop(2*ent1(ins1),R1(ins1),0,emax,bins);
    Dd  = histv0/Z;
    Ixx = histw0/Z;
    
    Rxy2 = U'*Rxy*V;
    R0 = -Rxy2 + Rxy2';
%    R1 = -Rxy2 - Rxy2';
    clearvars Rxy2 Rxy;
    % See Eq. 21 of 1507.01639
    R0 = pi*real( conj(R0).*R0 );
%    R1 = 2*pi*real( R1.*(R1.') );
    %R02 = bsxfun(@times,bsxfun(@times,R0,nf),nf');
    histw0 = ... 
        histwv_loop(reshape(ent0,S^2,1),reshape(R0,S^2,1),0,emax,bins,subs0,inds0);
    %R12 = bsxfun(@times,bsxfun(@times,R1,nf),1-nf');
%    [histw2, ~] = histwv_loop(2*ent1(ins1),R1(ins1),0,emax,bins,0);
    Ixy = histw0/Z;

    % The antisymmetric channel
    R0 = 1i*(-U'*RxyA*U + V'*RxyB*V);
    clearvars RxyA RxyB
%    R1 = 1i*(-R0 + R1);
    %= W0;
    % See Eq. 21 of 1507.01639
    R0 = pi*real( conj(R0).*R0 );
%    R1 = 2*pi*real( R1.*(R1.') );
    %R02 = bsxfun(@times,bsxfun(@times,R0,nf),nf');
    histw0 = ...
        histwv_loop(reshape(ent0,S^2,1),reshape(R0,S^2,1),0,emax,bins,subs0,inds0);
    %R12 = bsxfun(@times,bsxfun(@times,R1,nf),1-nf');
%    [histw2, ~] = histwv_loop(2*ent1(ins1),R1(ins1),0,emax,bins,0);
    Ixy2 = histw0/Z;

else

    %Assume T is a vector and make a matrix list of number operators for the 
        % filling of each quasiparticle
    %nf  = 1/( 1 + exp(bsxfun(@rdivide,en2,-T')) );

    %Store the filling fractions in a big matrix
    %nf0 = bsxfun(@times,repmat(nf,1,S),nf');
    %nf1 = bsxfun(@times,repmat(nf,1,S),1-nf');

    %Initialization
    %Dd = zeros(numel(T),bins); Ixx = Dd; Ixy = Dd; Ixy2 = Dd;
    
    ent1 = bsxfun(@minus,en2',en2);%-repmat(en2,1,S) + repmat(en2',S,1);
    %ins1 = (emax/2>=ent1)&(ent1>=0);

    %prepare fermi factors
%     nf0 = zeros(S,S,numel(T));
%     nf1 = zeros(S,S,numel(T));
%     for in = 1:numel(T)
%         if T(in) ~= 0
%             nf = 1./( 1 + exp(en2./(T(in))));
%             nf0(:,:,in) = bsxfun(@times,1-nf',1-nf);
%             nf1(:,:,in) = bsxfun(@times,1-nf',  nf);
%         else
%             nf0(:,:,in) = ones(S,S);
%             nf1(:,:,in) = ones(S,S);
%         end
%     end
    %reshape for convenience
    %nf0 = reshape(nf0,S^2,numel(T));
    %nf1 = reshape(nf1,S^2,numel(T));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Histogram it
    Rxx2 = U'*Rxx*V;
    R0 = -Rxx2 + Rxx2';
    R1 =  Rxx2 + Rxx2';
    clearvars Rxx2 Rxx
    % See Eq. 21 of 1507.01639
    R0 = single( abs(R0).^2 );
    R1 = single( 2*abs(R1).^2 );
    W1 = zeros(S^2,numel(T),'single');
    
    for in = 1:numel(T)
        if T(in) ~= 0
            nf = single( 1./( 1 + exp(en2./(T(in)))) );
            W0(:,in) = reshape( bsxfun(@times,bsxfun(@times,R0,1-nf'),1-nf) ,S^2,1);
            W1(:,in) = reshape( bsxfun(@times,bsxfun(@times,R1,1-nf'),nf) ,S^2,1);
        else
            W0(:,in) = reshape(R0,S^2,1);
        end
    end
    [histw0, histv0, subs0, inds0] = ...
        histwv_loop(reshape(ent0,S^2,1),W0,0,emax,bins);
    [histw1, histv1, subs1, inds1] = ...
        histwv_loop(reshape(ent1,S^2,1),W1,0,emax,bins);
    Dd  = (histv0+2*histv1)/Z;
    Ixx = (histw0+histw1)/Z;


    Rxy2 = U'*Rxy*V;
    R0 = -Rxy2 + Rxy2';
    R1 =  Rxy2 + Rxy2';
    clearvars Rxy2 Rxy;
    % See Eq. 21 of 1507.01639
    R0 = single( abs(R0).^2 );
    R1 = single( 2*abs(R1).^2 );
    for in = 1:numel(T)
        if T(in) ~= 0
            nf = single( 1./( 1 + exp(en2./(T(in)))) );
            W0(:,in) = reshape( bsxfun(@times,bsxfun(@times,R0,1-nf'),1-nf) ,S^2,1);
            W1(:,in) = reshape( bsxfun(@times,bsxfun(@times,R1,1-nf'),nf) ,S^2,1);
        else
            W0(:,in) = reshape(R0,S^2,1);
        end
    end
    histw0 = histwv_loop(reshape(ent0,S^2,1),W0,0,emax,bins,subs0,inds0);
    histw1 = histwv_loop(reshape(ent1,S^2,1),W1,0,emax,bins,subs1,inds1);
    Ixy = (histw0+histw1)/Z;

    % The antisymmetric channel
    R0 = U'*RxyA*U;
    R1 = V'*RxyB*V;
    temp = (R0 - R1);
    R1   = (R0 + R1);
    R0 = temp;
    clearvars RxyA RxyB
    % See Eq. 21 of 1507.01639
    R0 = single( abs(R0).^2 );
    R1 = single( 2*abs(R1).^2 );
    for in = 1:numel(T)
        if T(in) ~= 0
            nf = single( 1./( 1 + exp(en2./(T(in)))) );
            W0(:,in) = reshape( bsxfun(@times,bsxfun(@times,R0,1-nf'),1-nf) ,S^2,1);
            W1(:,in) = reshape( bsxfun(@times,bsxfun(@times,R1,1-nf'),nf) ,S^2,1);
        else
            W0(:,in) = reshape(R0,S^2,1);
        end
    end
    histw0 = histwv_loop(reshape(ent0,S^2,1),W0,0,emax,bins,subs0,inds0);
    histw1 = histwv_loop(reshape(ent1,S^2,1),W1,0,emax,bins,subs1,inds1);
    Ixy2 = (histw0+histw1)/Z;

end

%memory
%clearvars U V

%display the number of points above emax
%disp( sum( ceil( abs( en2(en2>emax/4) ) ) ) / ( size(en2,1)*size(en2,2) ) )
%disp( sum( ceil( abs( ent(ent>emax/2) ) ) ) / ( size(ent,1)*size(ent,2) ))

% I = cell(1,5);
% I{1} = Ev;    
% I{2} = Dd;
% I{3} = Ixx;   
% I{4} = Ixy;   
% I{5} = Ixy2; 

%toc
end
%  MIDAS SOLVER
%  Version 1.0.0 --- April 2016
%
%  Section of Biomedical Image Analysis
%  Department of Radiology
%  University of Pennsylvania
%  Richard Building
%  3700 Hamilton Walk, 7th Floor
%  Philadelphia, PA 19104
%
%  Web:   https://www.med.upenn.edu/sbia/
%  Email: sbia-software at uphs.upenn.edu
%
%  Copyright (c) 2018 University of Pennsylvania. All rights reserved.
%  See https://www.med.upenn.edu/sbia/software-agreement.html or COPYING file.

%  Author:
%  Erdem Varol
%  software@cbica.upenn.edu
%
%
% Reference: Varol, Erdem, Aristeidis Sotiras, Christos Davatzikos. 
% "MIDAS: regionally linear multivariate discriminative statistical mapping." NeuroImage (2018)



function map=midas_solver(X,Y,params,foreground)

if nargin<4
    foreground=ones(size(X{1}));
end
idx=find(foreground==1);
if nargin<3
    params={};
end
if ~isfield(params,'radius')
    params.radius=20;
end
if ~isfield(params,'num_neighborhoods')
    params.num_neighborhoods=100;
end
if ~isfield(params,'debug')
    params.debug=0;
end
if ~isfield(params,'mask_type');
    params.mask_type='spherical';
end
if ~isfield(params,'ls_svm_C');
    params.ls_svm_C=1;
end
if ~isfield(params,'ls_svm_type');
    params.ls_svm_type='erdem_slack_dual_woodbury';
end
%disp('MIDAS parameter settings:')
%disp(params)
template=X{1};
Xmat=zeros(length(X),numel(template)); %maybe save as sparse for memory?
for i=1:length(X)
    Xmat(i,:)=X{i}(:)';
end
n=length(X);
clear X
weights=ones(size(Y,1),1);
Y=zscore(Y); %Zero-mean unit-variance the labels
CY=cov(Y); %Covariance of labels
H=inv(CY); %inverse covariance of labels


%Pre allocating space for numerator and denominator of the MIDAS statistic for each covariate
S=cell(1,size(Y,2)); 
M=cell(1,size(Y,2));
for i=1:size(Y,2)
    S{i}=zeros(size(template));
    M{i}=zeros(size(template));
end

%Pre allocating space for cross covariance, variance and foreground coverage maps.
CC=cell(params.num_neighborhoods,1);
VV=zeros(size(template));
N=zeros(size(template));
N(foreground==0)=nan;

%Preallocating space for neighborhood masks
MASK0cell=cell(params.num_neighborhoods,1);
for t=1:params.num_neighborhoods
    if mod(t,1)==0
        disp(['iteration ' num2str(t)])
    end
    mask=mask_select(template,foreground,N,params); %Obtain local neighborhood mask
    mask=and(mask,foreground);
    [ii,~,vv]=find(mask(:));
    MASK0cell{t}=[ii repmat(t,length(ii),1) vv];
    maskData=bsxfun(@minus,Xmat(:,mask==1),mean(Xmat(:,mask==1),1)); %Z-score local data
    Q=dimensionality_reduction(maskData,100);
    Xq=maskData*Q; %Dimensionality reduction from d to n
    [w00,~,C00]=lssvm_3types_weighted(Xq,Y,weights,params); %Local LS-SVM solution
    w0=Q*w00;C0=Q*C00; %Dimensionality recovery from n to d
    C=(1/n)*maskData'*(maskData*C0); %Haufe correction to C matrix
    w=C*(Y*H); %Haufe correction to W to obtain activations
    for i=1:size(Y,2)
        S{i}(mask==1)=S{i}(mask==1)+w(:,i); %Additive update of the numerator of MIDAS statistic
        M{i}(mask==1)=M{i}(mask==1)+w0(:,i)'*w0(:,i); %Additive update of the denomerator of MIDAS statistic
    end
    CC{t}=sparse(repmat(ii,n,1),reshape(repmat((1:n),length(ii),1),n*length(ii),1),C(:),numel(template),n);
    N(mask==1)=N(mask==1)+1;
    VV(mask==1)=VV(mask==1)+sum(C0(:).^2);
end
MASK0IJV=cell2mat(MASK0cell);
clear MASK0cell
MASK0=sparse(MASK0IJV(:,1),MASK0IJV(:,2),MASK0IJV(:,3),numel(template),params.num_neighborhoods); %Store neighborhood mask locations in a sparse way to compute possible overlap later
clear MASK0IJV
MMT=MASK0'*MASK0;
clear MASK0


map.MMT=MMT;%Storing neighborhood overlap matrix
whos CC 
map.S=S;%Storing MIDAS numerator
map.M=M;%Storing MIDAS denominator
map.N=N;%Storing MIDAS coverage map
map.VV=VV;%Storing MIDAS variance
CV=zeros(numel(template),1);
disp('Covariance calculating');

%Below is how cross-covariance is computed
for p=1:params.num_neighborhoods
    disp([num2str(p) '/' num2str(params.num_neighborhoods)])
    for q=p:params.num_neighborhoods
        if MMT(p,q)%>MMT(1,1)*0.1
            tmp_pq=sum(CC{p}.*CC{q},2);
            if p~=q
                CV=CV+2*tmp_pq;
            else
                CV=CV+tmp_pq;
            end
        end
    end
end
map.CV=reshape(CV,size(template)); %storing MIDAS covarianace
P=cell(2,1);
for i=1:size(Y,2)
    map.stat{i}=(map.S{i}./map.M{i})./(sqrt((H(:,i)'*CY*H(:,i))*map.CV./map.VV.^2)); %Computing MIDAS stat
    P{i}=normcdf(-abs(map.stat{i}),0,1); %Computing p-value
    map.p{i}=reshape(2*P{i},size(template)); %Reshaping p-value to match image
end
map.params=params; %saving parameters
map.foreground=foreground; %saving foreground map
end

function coor=coordinate_select(template,foreground,N)
%Function to adaptively sample neighborhoods across the image
idx=find(and(foreground==1,N<=quantile(N(foreground==1),0.1))); %selecting neighborhoods that are undersampled
if isempty(idx)
    idx=find(foreground==1);
end
if size(template,3)>1
    r=randi(length(idx));
    [I,J,K]=ind2sub(size(template),idx(r));
    coor=[I J K];
elseif size(template,2)>1
    r=randi(length(idx));
    [I,J]=ind2sub(size(template),idx(r));
    coor=[I J];
end
end

function mask=mask_select(template,foreground,N,params)
%Function for obtaining a random neighborhood based on the desired topology
if strcmpi(params.mask_type,'spherical');
    coor=coordinate_select(template,foreground,N);
    mask=searchlight(template,coor,params.radius*(0.5+1*rand(1)));
elseif strcmpi(params.mask_type,'random');
    coor=coordinate_select(template,foreground,N);
    mask=searchlight(template,coor,params.radius*(0.5+1*rand(1)));
    mask(randperm(numel(mask)))=mask;
elseif strcmpi(params.mask_type,'bispherical');
    coor1=coordinate_select(template,foreground,N);
    coor2=coordinate_select(template,foreground,N);
    mask1=searchlight(template,coor1,(1/sqrt(2))*params.radius*(0.5+1*rand(1)));
    mask2=searchlight(template,coor2,(1/sqrt(2))*params.radius*(0.5+1*rand(1)));
    mask=or(mask1,mask2);
elseif strcmpi(params.mask_type,'trispherical');
    coor1=coordinate_select(template,foreground,N);
    coor2=coordinate_select(template,foreground,N);
    coor3=coordinate_select(template,foreground,N);
    mask1=searchlight(template,coor1,(1/sqrt(3))*params.radius*(0.5+1*rand(1)));
    mask2=searchlight(template,coor2,(1/sqrt(3))*params.radius*(0.5+1*rand(1)));
    mask3=searchlight(template,coor3,(1/sqrt(3))*params.radius*(0.5+1*rand(1)));
    mask=or(or(mask1,mask2),mask3);
elseif strcmpi(params.mask_type,'quadspherical');
    coor1=coordinate_select(template,foreground,N);
    coor2=coordinate_select(template,foreground,N);
    coor3=coordinate_select(template,foreground,N);
    coor4=coordinate_select(template,foreground,N);
    mask1=searchlight(template,coor1,(1/sqrt(4))*params.radius*(0.5+1*rand(1)));
    mask2=searchlight(template,coor2,(1/sqrt(4))*params.radius*(0.5+1*rand(1)));
    mask3=searchlight(template,coor3,(1/sqrt(4))*params.radius*(0.5+1*rand(1)));
    mask4=searchlight(template,coor4,(1/sqrt(4))*params.radius*(0.5+1*rand(1)));
    mask=or(or(mask1,mask2),or(mask3,mask4));
end
end

function [w,b,C]=lssvm_3types_weighted(X,Y,S,params)
%Least squares SVM solver along with hat matrix C
type=params.ls_svm_type;
c=params.ls_svm_C;
idx=find(S>0);
[n0,~]=size(X);
X=X(idx,:);
Y=Y(idx,:);
S=S(idx,:);
[n,d]=size(X);
c=c*d/n;
%% bilwaj
if strcmpi(type,'bilwaj')
    warning('ls svm type not supported')
    iXXt=inv(X*X');
    J=ones(n,1);
    C=X'*(iXXt+iXXt*J*inv(-J'*iXXt*J)*J'*iXXt);
    w=C*Y;
    b=-mean(X*w-Y,1);
end
%% me

if strcmpi(type,'erdem')
    warning('ls svm type not supported')
    A=[eye(d) zeros(d,1) X'; zeros(1,d) zeros(1,1) ones(1,n); X ones(n,1) zeros(n,n)];
    B=[zeros(d,1);zeros(1,1);Y];
    % v=inv(A'*A)*A'*B;
    C0=inv(A'*A)*A';
    v=C0*B;
    C=C0(1:d,d+2:end);
    w=v(1:d,1);
    b=v(d+1);
end

%% me with slack KKT

if strcmpi(type,'erdem_slack')
    A=[eye(d) zeros(d,1) zeros(d,n) X';zeros(1,d) zeros(1,1) zeros(1,n) ones(1,n); zeros(n,d) zeros(n,1) -eye(n) c*diag(S);X ones(n,1) -eye(n) zeros(n,n)];
    B=[zeros(d,1);zeros(1,1);zeros(n,1);Y];
    v=A\B;
    w=v(1:d,1);
    b=v(d+1);
    C0=inv(A);
    C=C0(1:d,d+1+n+1:end);
end

%% me with slack + dual

if strcmpi(type,'erdem_slack_dual')
    A=[zeros(1,1) ones(1,n);ones(n,1) -X*X'-inv(diag(S))/c];
    B=[zeros(1,1);Y];
    v=A\B;
    lambda=v(2:end);
    w=-X'*lambda;
    b=v(1);
    C0=inv(A);
    C=-X'*C0(2:end,2:end);
end

%% me with slack + dual + woodbury inversion

if strcmpi(type,'erdem_slack_dual_woodbury')
    K=inv(X*X'+inv(diag(S))/c);
    J=ones(n,1);
    C=X'*(K+(K*J/(-J'*K*J))*J'*K);
    w=C*Y;
    b=J'*(Y-X*w)/n;
end

%% me with slack + primal

if strcmpi(type,'erdem_slack_primal')
    A=[eye(d)/c+X'*diag(S)*X X'*diag(S)*ones(n,1);ones(1,n)*diag(S)*X ones(1,n)*diag(S)*ones(n,1)];
    B=[X'*diag(S);ones(1,n)*diag(S)]*Y;
    v=A\B;
    w=v(1:d);
    b=v(d+1);
    C0=inv(A)*[X';ones(1,n)];
    C=C0(1:d,:);
end

%% Suykens

if strcmpi(type,'suykens')
    error('ls svm type not supported')
    Z=bsxfun(@times,X,Y);
    A=[eye(d) zeros(d,1) zeros(d,n) -Z'; zeros(1,d) zeros(1,1) zeros(1,n) -Y';...
        zeros(n,d) zeros(n,1) c*eye(n) -eye(n);Z Y eye(n) zeros(n)];
    B=[zeros(d,1);zeros(1,1);zeros(n,1);ones(n,1)];
    v=inv(A)*B;
    w=v(1:d,1);
    b=v(d+1);
end

CC=zeros(d,n0);
CC(:,idx)=C;
C=CC;clear CC
end

function mask=searchlight(X,coor,radius)
%Function for generating a neighborhood mask based on the coordinates and radius given, supports 2D or 3D images
if size(X,3)>1
    [x,y,z] = meshgrid(1:size(X,2),1:size(X,1),1:size(X,3));
    mask = sqrt((x-coor(2)).^2 + (y-coor(1)).^2 + (z-coor(3)).^2)<=radius;
elseif size(X,2)>1
    [x,y]=meshgrid(1:size(X,2),1:size(X,1));
    mask = sqrt((x-coor(2)).^2 + (y-coor(1)).^2)<=radius;
end
end

function Q=dimensionality_reduction(X,dim)
%Function for dimensionality reduction based on random projections
%Input: X: N x D data matrix
%Input: dim: reduced dimensionality
%Output: Xreduced: N x dim approximate data matrix
%Output: Q: reconstruction orthonormal matrix
[N,~]=size(X);
omega=randn(N,dim);
omega=bsxfun(@times,omega,1./norms(omega,2,1));
Y=X'*omega;
clear omega
% [Q,~]=mygramschmidt(Y);
[Q,~]=qr(Y,0);
clear Y
%Xreduced=(Q'*X')';
end

function o = norms( x, p, dim )
% Function to compute vector norms
switch p,
    case 1,
        o = sum( abs( x ), dim );
    case 2,
        o = sqrt( sum( x .* conj( x ), dim ) );
    case Inf,
        o = max( abs( x ), [], dim );
    otherwise,
        o = sum( abs( x ) .^ p, dim ) .^ ( 1 / p );
end

end

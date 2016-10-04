function [rank_scr,rank_cho,rank_scr_supl]=rank_scr_cho_supl(A,label)
%   Identifying differentially expressed genes by scotter_compact_ratio ----> scr, Using formula_2: delta_2 and d_2 ( in paper) ,which euqalt to delta_3 and d_3  in funcion  'yk_DS_GeneFilter.m'
%   Identifying differentially expressed genes by function of FEBS_cho ----> rank_cho et al
%   Identifying differentially expressed genes by scotter_compact_ratio ----> d_1+delta_1 (in paper)
%   
%   
%  -------------------------------------------Version 1.0----------------------------------------
%       Usage: rank_scr_cho_supl(train_data,train_label)                 
%       Input Parameters:   A       ------ training samples    row is samples ; column is genes
%                           label   ------ training smaples'labels, positive integer, is column vector
%       Output: 
%           rank_scr is a matrix in size of number_of_gene*2, the first column is the Index of gene, the second column is sorted scotter_comp_ratio
%           rand_cho is a matrix in size of number_of_gene*2, the first column is the Index of gene, the second column is sorted cho_score
%           rank_scr_supl is a matrix in size of number_of_gene*2, the first column is the Index of gene, the second column is also another sorted scotter_comp_ratio
%       yang kun programmed  2005-9-3
%  -----------------------------------------------------------------------------------------------

epsilon=1e-6; % 2006-1-19 % Prevent some elements of scatter equalling to 0;

if (nargin <=1) 
    fprintf(1,'\n  Input Error: no training data or no training data label.\n\n');
    help rank_scr_cho_supl.m
    return
end
[nTotalSample,nTotalGene]=size(A);
nTotalClass=max(label);

A_centroid=zeros(nTotalClass,nTotalGene);   % centroid of each class about data matrix
A_centroid_mean=zeros(1,nTotalGene);    % mean of centroids

X=zeros(nTotalSample,nTotalGene);   % difference matrix
X_centroid=zeros(nTotalClass,nTotalGene);   % centroid of each class about difference matrix
X_centroid_mean=zeros(1,nTotalGene);    % mean of all class centroids

for i=1:nTotalClass
    index=find(label==i);
    A_centroid(i,:)=mean(A(index,:),1);
end

% circle is made by vector
for i=1:nTotalClass
    index=find(label==i);
    % duplicate row vector "A_centroid(i,:)" into length(index) rows, size=[length(index),nTotalGene] ;
    X(index,:)=abs(A(index,:)-A_centroid(i*ones(length(index),1),:)); 
    X_centroid(i,:)=mean(X(index,:),1);
end

A_centroid_mean=mean(A_centroid,1);
X_centroid_mean=mean(X_centroid,1);

% ---------computing the variance "S","scatter", which are equality in "formula_2" and "formula_3"--------------
S=zeros(1,nTotalGene);
scatter=zeros(1,nTotalGene);

S=std(A_centroid,1,1);
if (nTotalClass==2)
    scatter=S+abs(A_centroid(1,:)-A_centroid(2,:));
elseif(nTotalClass>2)
    A_centroid=sort(A_centroid,1);          %  Attention !!! After here, "A_centroid" is differnet with normal 
    scatter=S+0.5*min(A_centroid(2:nTotalClass,:)-A_centroid(1:nTotalClass-1,:),[],1);
end


% --------------------- computing the variance "delta_2","d_2","compact_2","score_2", in "formula_2" ----------------
d_2=zeros(1,nTotalGene);
delta_2=zeros(1,nTotalGene);
compact_2=zeros(1,nTotalGene);
score_2=zeros(1,nTotalGene);

for i=1:nTotalClass
    index=find(label==i);
    d_2=sum(X(index,:).^2,1)/length(index)+d_2;
end
d_2=d_2/nTotalClass;
delta_2=d_2-X_centroid_mean.^2;

d_2=sqrt(d_2);
delta_2=sqrt(delta_2);

compact_2=delta_2+d_2;

    %%% 2006-1-19 revised
    index_special=find(scatter==0);
    scatter(index_special)=epsilon; % Prevent some elements of scatter equalling to 0;
    %%% 2005-1-19 revised
% computing the variance "score_2"
score_2=compact_2./scatter;

rank_scr=zeros(2,nTotalGene);
rank_scr(1,:)=1:nTotalGene;
rank_scr(2,:)=score_2;
rank_scr=rank_scr';
rank_scr=sortrows(rank_scr,2);


% --------------------- computing the variance "delta_1","d_1","compact_1","score_1", in "formula_1"
d_1=zeros(1,nTotalGene);
delta_1=zeros(1,nTotalGene);
compact_1=zeros(1,nTotalGene);
score_1=zeros(1,nTotalGene);

d_1=sum(X_centroid.^2,1)/nTotalClass;
d_1=sqrt(d_1);
delta_1=std(X_centroid,1,1);

compact_1=delta_1+d_1;

% computing the variance "score_1"
score_1=compact_1./scatter;

rank_scr_supl=zeros(2,nTotalGene);
rank_scr_supl(1,:)=1:nTotalGene;
rank_scr_supl(2,:)=score_1;
rank_scr_supl=rank_scr_supl';
rank_scr_supl=sortrows(rank_scr_supl,2);


% --------------------------- computing the "FEBS_gene_score" -------------------------------------------------
FEBS_compact=zeros(1,nTotalGene);
for i=1:nTotalClass
    index=find(label==i);
    X(index,:)=X(index,:)-X_centroid(i*ones(length(index),1),:);
end
FEBS_compact=sum(X.^2,1);
FEBS_compact=FEBS_compact.*nTotalSample;
FEBS_compact=FEBS_compact./((nTotalSample-1)*nTotalClass);
FEBS_compact=sqrt(FEBS_compact);
FEBS_compact=FEBS_compact.*X_centroid_mean;

    %%% 2006-1-19 revised
    index_special=find(S==0);
    S(index_special)=epsilon; % Prevent some elements of scatter equalling to 0;
    %%% 2005-1-19 revised

FEBS_score=FEBS_compact./S;

rank_cho=zeros(2,nTotalGene);
rank_cho(1,:)=1:nTotalGene;
rank_cho(2,:)=FEBS_score;
rank_cho=rank_cho';
rank_cho=sortrows(rank_cho,2);
return
% function is OVER

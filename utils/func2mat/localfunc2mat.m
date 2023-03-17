function [A,ysiz]=localfunc2mat(fcn, xsiz,widths)
%LOCALFUNC2MAT - Efficient algorithm for obtaining the matrix representing a 
%linear function from one N-D image space to another, given that each pixel
%of the input image has a very localized effect on the output image.
%
%
%TERMINOLOGY: We are given a linear transformation Y=f(X), where X is an
%N-dimensional input image, and Y is an M-dimensional output image. Suppose
%that two input images X1 and X2 have zero elements except for a single element
%in each, X1(i)=1 and X2(j)=1. Then we say that w is the "influence width"
%of f() if sufficient separation in i,j that is |i-j|>w, implies that
%Y1=f(X1) and Y2=f(X2) have no overlap in their non-zero regions. Here i,j
%are vector-valued coordinates and w is either a scalar or a vector. If w
%is vector-valued, the  influence width is different along different
%dimensions of X. In the case where the output image space is the same as the 
%input image space, the influence widths are the same as the point spread 
%function (PSF) widths.
%
%
%CODE SYNTAX: The input syntax is,
%
%  [A,ysiz]=localfunc2mat(fcn, xsiz,inwidths)                           
% 
%where,
%
%IN:                         
%                            
%       fcn: A function handle to f() accepting the syntax fcn(X).                
%      xsiz: The known input image dimensions, size(X).                
%    inwidths: The vector of influence widths. Can also be a scalar if
%              the influence widths are the same for all dimensions
%             
%
%OUT:                        
%                            
%       A:  The output sparse matrix. The dimensions of A are MxN where 
%           N=numel(X) and M=numel(fcn(X)).
%
%    ysiz:  The output dimensions, size(fun(X))


X=zeros(xsiz,'like',xsiz);
xsiz=size(X);

if isscalar(widths)
   widths(1:numel(xsiz))=widths;
else
    assert(numel(widths)==numel(xsiz),'WIDTHS must be same size as SIZ or a scalar')
end

subsets=widths+2;

N=numel(subsets);

starts=arrayfun(@(c) 1:c,subsets,'uni',0);

[C{1:N}]=ndgrid(starts{:});

combinations =reshape(  cat(N+1,C{:})   , [],N  ); clear C

K=size(combinations,1);

Cols=reshape(1:prod(xsiz), xsiz);

[Icell, Jcell, Scell]=deal(cell(K,1));

for k=1:K

    idxChequer=arrayfun(@(a,b,c) a:b:c, combinations(k,:),subsets,xsiz,'uni',0);
    
      X(idxChequer{:})=1;
    
     Y=fcn(X);

      X(idxChequer{:})= Cols(idxChequer{:});

    L=bwlabeln(Y);
    
    I=find(L~=0);
           
    G=L(I);
    S=Y(I);
    clear Y
    
    denom=splitapply(@sum, S.^2, G);
    
          
      Yw=fcn(X);    
    num=splitapply(@sum, S.*Yw(I),  G);
      clear Yw
    
    j=round(num./denom);
    J=j(G);
    
    [Icell{k}, Jcell{k}, Scell{k}]=deal(I,J,S);
    
    X(:)=0;
end
    
  I=cell2mat(Icell); J=cell2mat(Jcell); S=cell2mat(Scell);
  
  M=numel(L);
  N=numel(X);
  ysiz=size(L);
  
  A=sparse(I,J,S,M,N);
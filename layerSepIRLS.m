function [I1,I2] = layerSepIRLS (I,G, f_inds,b_inds)
%Layers reconstruction using Iterative Reweighted Least Square (IRLS)
%   [ I1 I2 ] = layerSepIRLS( I, G, edgeProbB, edgeProbF)
%   I: mixture image
%   G: derivative filter matrix (pre-computed)
%   f_inds,b_inds:  pixels labelled as belonging to reflection/background

iter = 5;

w1=50;
w2=1;

[h,w]=size(I);
imgSize=h*w;

%constraints matrix
Gx = G.Gx;
Gy = G.Gy;
Gxx = G.Gxx;
Gyy = G.Gyy;


A1=[Gx;Gy;Gxx;Gyy];
A=[A1;A1];
b=[zeros(size(A1,1),1);A1*I(:)];
f1=[ones(imgSize*2,1)*w1;w2*ones(imgSize*2,1);w2*ones(imgSize,1)];

f2=f1;
f1([f_inds,f_inds+imgSize])=0;
f1([b_inds,b_inds+imgSize])=100;
f2([f_inds,f_inds+imgSize])=100;
f2([b_inds,b_inds+imgSize])=0;

f1([b_inds+imgSize*2,b_inds+imgSize*3,b_inds+imgSize*4])=4;
f2([f_inds+imgSize*2,f_inds+imgSize*3,f_inds+imgSize*4])=4;

f=[f1;f2];

rinds1=find(sum(A~=0,2)==0);
rinds2=find(f==0);
inds=setdiff([1:size(A,1)],[rinds1;rinds2]);
A=A(inds,:);
b=b(inds,:);
f=f(inds);
oA=A;
ob=b;

df=spdiags(f,0,length(f),length(f));
A=df*oA; b=df*ob;
x=(A'*A)\(A'*b);

fprintf('Initial error = %g \n',sum(abs(A*x-b)));
for j=1:iter
      e=abs(oA*x-ob);
      e=max(e,0.00001);
      e=1./e;
      E=spdiags(e,0,length(f),length(f));
      x=(A'*E*A)\(A'*E*b);
      fprintf('iteration = %d, current residue = %g \n', j, sum(abs(A*x-b)));   
end

I1=reshape(x,h,w);
I2=I-I1;



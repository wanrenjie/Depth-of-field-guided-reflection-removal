clear all
img = imread('153.jpg');
cform = makecform('srgb2lab');
I = applycform(img,cform);

L = I(:,:,1); a = I(:,:,2); b = I(:,:,3);
% % 
PyramidL = BuildPyramid(L);
Pyramida = BuildPyramid(a);
Pyramidb = BuildPyramid(b);

%Fusion
% a channel
[N1,N2,N3] = size(PyramidL{1});
pl1 = PyramidL{1};
pl2 = imresize(PyramidL{2},[N1,N2]);
pl3 = imresize(PyramidL{3},[N1,N2]);

lmap = (0.4*pl2+0.6*pl3).*(pl1);
Thl = Threshold(lmap);
map1 = heaviside(lmap - Thl);

pa1 = Pyramida{1};
pa2 = imresize(Pyramida{2},[N1,N2]);
pa3 = imresize(Pyramida{3},[N1,N2]);

amap = (0.4*pa2+0.6*pa3).*(pa1);
Tha = Thl/1.5;
map2 = heaviside(amap - Tha);

%b channel
pb1 = Pyramidb{1};
pb2 = imresize(Pyramidb{2},[N1,N2]);
pb3 = imresize(Pyramidb{3},[N1,N2]);

bmap = (0.4*pb2+0.6*pb3).*(pb1);
Thb = Tha/1.5;
map3 = heaviside(bmap - Thb);

Background = map1|map2|map3;

%%obtain the reflection component
se = strel('line',10,10);

w = fspecial('sobel');
for ch = 1:size(img,3)
    img = im2double(img);
    gx(:,:,ch) = imfilter(img(:,:,ch),w);
    gy(:,:,ch) = imfilter(img(:,:,ch),w');
end
grad = sqrt(gx.^2+gy.^2);
grad = max(grad,[],3);
reflectionPoints = find(grad<0.3&grad>0.05);
map4 = zeros(N1,N2);
map4(reflectionPoints) = 1;

bw2 = map4;
reflectionpoints = find(bw2 == 1);
for i = 1:length(reflectionpoints)
    p = reflectionpoints(i);
    [Hp,rows,cols] = getpatch([N1,N2],p);
    o = Background(Hp);
    flag = find(o==1);
    if(Background(p)==1|length(flag)~=0)
        bw2(p) = 0;
    end
end
% Reflection = imdilate(bw2,se);
Reflection = bw2;
indF = find(Reflection == 1);
indB = find(Background == 1);

[h,w,d] = size(grad);
G = struct;
[G.Gx,G.Gy,G.Gxx,G.Gxy,G.Gyy]=getGMat(w,h);

for ch = 1:3
    [I1(:,:,ch),I2(:,:,ch)]=layerSepIRLS (img(:,:,ch),G,indF,indB);
end
figure
imshow(I1);
figure
imshow(I2);






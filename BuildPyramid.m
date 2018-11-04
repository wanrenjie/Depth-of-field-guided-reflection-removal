function [ Pyramid ] = BuildPyramid(img)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%     cform = makecform('srgb2lab');
    nbins = 41;
    k = [3,5,7];
%     I = applycform(img,cform);
%     L = I(:,:,1); a = I(:,:,2);b = I(:,:,3);

    L = im2double(img);
    f1 = [1,-1];f2 = [1;-1];

    rp = [1,0.8,0.5];%Resize Parameter
    Pyramid = cell(1,3);
%   DoF = cell(1,3);
    D = cell(1,3);
    for kk = 1:3
        im = imresize(L,rp(kk));
        
        rhox1 = -imfilter(im,f1,'circular');
        rhox1 = (rhox1+1)./2;
        px1 = makeHistogram(nbins,rhox1);
    
        rhoy1 = -imfilter(im,f2,'circular');
        rhoy1 = (rhoy1+1)./2;
        py1 = makeHistogram(nbins,rhoy1);
    
        for ii = 1:3
            [N1,N2] = size(im);
            G = fspecial('gaussian',[k(ii) k(ii)],5);       
            rhoxk = imfilter(im,G,'same');
            rhoyk = imfilter(im,G,'same');
        
            rhoxk =  -imfilter(rhoxk,f1,'circular');
            rhoxk = (rhoxk+1)./2; %This is very important, Otherwise we cannot get what we want.
            pxk = makeHistogram(nbins,rhoxk);
        
            rhoyk =  -imfilter(rhoyk,f2,'circular');
            rhoyk = (rhoyk+1)./2; %This is very important, Otherwise we cannot get what we want.
            pyk = makeHistogram(nbins,rhoxk);

            map = zeros(N1,N2);
            for y = 2 : N1-1
                for x = 2 : N2-1
                    map(y,x) = CalculateLogLikehood(x,y,rhoxk,rhox1,rhoyk,rhoy1,pxk,px1,pyk,py1,nbins);
                end
            end
            D{ii} = map/70; 
        end
        DoF = D{1}+D{2}+D{3};
        Pyramid{kk} = DoF;
    end
    
    function h = makeHistogram(nbins,img)   
        sum = 0;
        h = zeros(1,nbins);
        [height,width] = size(img);
        for yy = 1:height
            for xx = 1:width
                v = img(yy,xx);
                bin = uint8(v*nbins);
                if(bin>=nbins)
                    bin = bin-1;
                end
                h(bin) = h(bin)+1;
                sum = sum+1;
            end
            
        end
        h(find(h==0)) = 0.00000001;
        h = h./sum;

    function LL = CalculateLogLikehood(x,y,dxk,dx1,dyk,dy1,pxk,px1,pyk,py1,nbins)
        border = 2;
        LL = 0;
        for ii = y-border+1:y+border-1
            for jj = x-border+1:x+border-1
                vxk = dxk(ii,jj);
                bxk = uint8(vxk*nbins);
                if(bxk >= nbins) 
                    bxk = bxk - 1;
                end
                
                vx1 = dx1(ii,jj);
                bx1 = uint8(vx1*nbins);
                if(bx1 >= nbins)
                    bx1 = bx1 - 1;
                end
                
                vyk = dyk(ii,jj);
                byk = uint8(vyk*nbins);
                if(byk >= nbins) 
                    byk = byk - 1;
                end
                
                vy1 = dy1(ii,jj);
                by1 = uint8(vy1*nbins);
                if(by1 >= nbins)
                    by1 = by1 - 1;
                end
                
                pxk(bxk) = pxk(bxk)/(pxk(bxk) + px1(bx1));
                px1(bx1) = px1(bx1)/(pxk(bxk) + px1(bx1));
                
                pyk(byk) = pyk(byk)/(pyk(byk) + py1(by1));
                py1(by1) = py1(by1)/(pyk(byk) + py1(by1));
                
                
                LL = LL + (pxk(bxk)*log(pxk(bxk)/px1(bx1)) + pyk(byk)*log(pyk(byk)/py1(by1)));
%                 LL = LL + (pxk(bxk)*(pxk(bxk)/px1(bx1)) + pyk(byk)*(pyk(byk)/py1(by1)));
            end
        end


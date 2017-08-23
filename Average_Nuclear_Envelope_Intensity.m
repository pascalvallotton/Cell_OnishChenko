function [IntNE,IntNEMax,IntNEBack,IntNEBackInside,IntNEBackOutside,postArrayArrayArray,fileList] = Average_Nuclear_Envelope_Intensity()
%
% Average_Nuclear_Envelope_Intensity measures the
% average intensity (on an image- by image basis) in the Nuclear envelope. Intensities are
% measured on background-subtracted images. The transmittance image is used to exclude area that are out of focus.
% The intensity immediately around the detected nuclear envelope is used to
% precisely estimate ring intensities. Ring intensities are measured as
% averages of Nups's channel gated from the DHEL channel. Radius of NUP
% ring is estimated relative to HDEL ring.
%

% questions -> Pascal, Sasi & Evgeny

% Output IntNE contains all average intensities for ring for each image of the
% sequence. IntNEMax retains the maximum rim intensities. IntNEBack provides an estimate of the signal in the
% immediate neighbourhood of the mask (free nups)

% parameters
YeastMeanDiam=100; %Modify only when using different microscope... in pixels

pathName=uigetdir();drawnow;

selectedImage = rdir([pathName,'\**\*.nd2']);

postArrayArrayArray=[];

for numIm=1:size(selectedImage,1) %loop over all selected images
    try
        r = bfopen([selectedImage(numIm).name]);
        if strfind(selectedImage(numIm).name,'Mask') %don't analyse mask image
            continue
        end
        
        if ~iscell(selectedImage)
            ImWide=r{1, 1}{1, 1};
            ImFlu2=r{1, 1}{2, 1};
            ImFlu1=r{1, 1}{3, 1};
        else %load single image
            ImWide=Gauss2D(double(bfimread([pathName,selectedImage{numIm}],1)),1);
            ImFlu2=Gauss2D(double(bfimread([pathName,selectedImage{numIm}],2)),1);%figure;imshow(ImWide,[])
            ImFlu1=Gauss2D(double(bfimread([pathName,selectedImage{numIm}],3)),1);
        end
        ImWide=double(ImWide);ImFlu1=double(ImFlu1);ImFlu2=double(ImFlu2);
        
        % Identify out of focus area in widefield to exclude dodgy cells...You may have to increase the quantile if
        % you want to be more selective regarding removing area containing out
        % of focus cells
        threshQuantile = quantile(ImWide(:),0.999);
        
        % Create larger mask around blurred area to exclude questionable zones
        % (overlaying, out of focus etc) cells .
        
        maskExcl=imdilate(bwareaopen((ImWide>threshQuantile),(YeastMeanDiam/20)^2),strel('disk',round(YeastMeanDiam/2)));%figure;imshow(ImWide,[min(ImWide(:)),max(ImWide(:)/3)]);hold on; contour(maskExcl,'r','LineWidth',2)
        
        % background subtraction. ...removes bias caused by shading etc
        % we are processing at the same time both the HDEL channel (1...) and the NUP channel (2...)
        
        backFlu1=Gauss2D(imtophat(ImFlu1,strel('disk',round(YeastMeanDiam/10))),0.5);
        backFlu2=Gauss2D(imtophat(ImFlu2,strel('disk',round(YeastMeanDiam/10))),0.5);%figure;imshow(backFlu1,[])
        
        threshQuantile = quantile(backFlu1(:),0.960); % threshQuantile is used to select only bright bloobs among all blobs yet to be created
        threshQuantile2 = quantile(backFlu2(:),0.980);
        
        % maskNE is created in order to let the image outer boundary propagate to the
        % cell cores via an open path. ... otherwise most cells would be missed
        maskNE=imdilate(bwareaopen((backFlu1>threshQuantile),2),strel('disk',round(YeastMeanDiam/5)));% figure,imshow(maskNE,[]);figure;imshow(backFlu1>threshQuantile,[])
        maskNE2=imdilate(bwareaopen((backFlu2>threshQuantile2),2),strel('disk',round(YeastMeanDiam/5)));% figure,imshow(maskNE2,[]);figure;imshow(backFlu1,[])
        
        % Find all candidate regions delimited by edges
        edgeplot=edge((backFlu1),'log',0).*(1-maskExcl).*maskNE;
        edgeplot2=edge((backFlu2),'log',0).*maskNE2;
        
        edgeplotCut=bwmorph(edgeplot,'bridge');% cut into smaller pieces so you do not lose pieces due to region averaging (a dark region connected by a thin gap to a bright region obliterates entire region)
        labelMask=bwlabeln(1-edgeplotCut,4);
        regPropMask=regionprops(labelMask,backFlu1,'Area','MeanIntensity');
        idx = find([regPropMask.MeanIntensity] > 0.9*threshQuantile); %dim HDEL rings are removed
        idxExl = find([regPropMask.Area] > 10000); %giant bloobs are removed
        BW2 = ismember(labelMask, idx);%figure;imshow(4009*BW2+backFlu1,[])
        BW2Excl = ismember(labelMask, idxExl);
        
        edgeplotCut2=bwmorph(edgeplot2,'bridge');% cut into smaller pieces so you do not lose pieces due to region averaging
        labelMask2=bwlabeln(1-edgeplotCut2,4);
        regPropMask2=regionprops(labelMask2,backFlu2,'Area','MeanIntensity');
        idx2 = find([regPropMask2.MeanIntensity] > threshQuantile2);
        idxEx2 = find([regPropMask2.Area] > 10000);
        BW22 = ismember(labelMask2, idx2);
        BW2Excl2 = ismember(labelMask2, idxEx2);%figure;imshow(10900*BW22+backFlu1,[])
        
        BW2=bwareaopen(BW2,8); %remove small bits
        BW22=bwareaopen(BW22,8);
        
        labelMask=bwlabeln(1-BW2,8);
        labelMask2=bwlabeln(1-BW22,8);
        
        labelMask(find(BW2Excl))=0;
        regPropMask=regionprops(labelMask);
        idx = find([regPropMask.Area] < 8);
        BWx = ismember(labelMask, idx);
        BW2=BW2+BWx;
        
        labelMask2(find(BW2Excl2))=0;
        regPropMask2=regionprops(labelMask2);
        idx2 = find([regPropMask2.Area] < 8);
        BWx2 = ismember(labelMask2, idx2);
        BW22=BW22+BWx2;
        
        %remove ring components too wide to correspond to a bona fide ring
        TestWidthIm=bwdist(1-BW2); %distance transform
        labelMask=bwlabeln(TestWidthIm,4);
        regPropTest=regionprops(labelMask,TestWidthIm,'MaxIntensity');
        idxTest = find([regPropTest.MaxIntensity] >= 5.1);
        BW2T = ismember(labelMask, idxTest); %figure;imshow(10000*BW2T+ImFlu1,[])
        BW2=BW2.*(1-BW2T);%figure;imshow(-2000*BW2+backFlu1,[])
        
        TestWidthIm2=bwdist(1-BW22);
        labelMask2=bwlabeln(TestWidthIm2,4);
        regPropTest2=regionprops(labelMask2,TestWidthIm2,'MaxIntensity');
        idxTest2 = find([regPropTest2.MaxIntensity] >= 5.1);
        BW2T2 = ismember(labelMask2, idxTest2);
        BW22=BW22.*(1-BW2T2);
        
        % HDEL rings are sometimes broken you need to restore continuity. First locate
        % entire cell using ER adjacent to cell wall
        BW4=imclose(BW2,strel('disk',5));
        BW5=imfill(BW4,'holes');
        BW5=((BW5-BW4)>0);
        labelMask=bwlabeln(BW5,4);
        
        %Do the same with the NUPs
        BW42=imclose(BW22,strel('disk',5));
        BW52=imfill(BW42,'holes');
        BW52=((BW52-BW42)>0);
        labelMask2=bwlabeln(BW52,4);%blinkOrigResult( BW5,backFlu1,30,0,1)
        
        regPropMask=regionprops(labelMask);%figure;imshow(-2000*BW52+backFlu2,[])
        idx = find(or([regPropMask.Area] > 1100,[regPropMask.Area] < 270));
        BWx = ismember(labelMask, idx);
        BW7=(BW5-BWx>0);
        
        regPropMask2=regionprops(labelMask2);
        idx2 = find(or([regPropMask2.Area] > 1100,[regPropMask2.Area] < 270));
        BWx2 = ismember(labelMask2, idx2);
        BW72=(BW52-BWx2>0);
        
        BW8=imdilate(BW7,strel('disk',5))-BW7;
        BW82=imdilate(BW72,strel('disk',5))-BW72;
        labelMask=bwlabeln(BW7,8);
        thinNE=zeros.*BW4;
        regPropMask=regionprops(labelMask,backFlu2,'MeanIntensity','BoundingBox','Image','Centroid');
        postArrayArray=[];
        
        for k=1:size(regPropMask,1)
            try
                a=regPropMask(k).BoundingBox(2)-8;
                b=regPropMask(k).BoundingBox(1)-8;
                c=regPropMask(k).BoundingBox(4)+16;
                d=regPropMask(k).BoundingBox(3)+16;
                subimg=backFlu1(a:a+c,b:b+d);
                subimgFlu2=backFlu2(a:a+c,b:b+d);
                subimgFlu2Orig=ImFlu2(a:a+c,b:b+d);
                
                subimgBW8=BW8(a:a+c,b:b+d);
                subimg2=subimg.*0;
                subimg2(9:8+regPropMask(k).BoundingBox(4),9:8+regPropMask(k).BoundingBox(3))=regPropMask(k).Image;
                subimg3=imdilate(subimg2,strel('disk',4))-imdilate(subimg2,strel('disk',2));
                meanInt=mean(subimgFlu2(find(subimg3)));IntegratedNPCsignal=sum(subimgFlu2(find(subimg3)));
                posMax=LocMaximas1D(subimgFlu2(find(subimg3)),5);
                intsVals=subimgFlu2(find(subimg3));
                meanIntLocMax=quantile((intsVals(posMax)),0.8);
                subimg4=imdilate(subimg3,strel('disk',8))-imdilate(subimg3,strel('disk',6));
                
                subimg4Inside=(subimg4&subimg2);
                subimg4Outside=(subimg4-bwmorph(subimg2,'dilate'))>0;
                
                thinM=subimg3;
                meanIntBack=mean(subimgFlu2(find(subimg4)));
                meanIntBackInside=mean(subimgFlu2Orig(find(subimg4Inside)));
                meanIntBackOutside=mean(subimgFlu2Orig(find(subimg4Outside)));
                %does it overlap with a neighbour? if so ignore the data point to be on the safe side
                if max(max(subimgBW8+subimg4))>1
                    continue
                end
                %does it have any significant signal from nups? if nothing at all is present, we do not consider the data point, in order to be on
                %safe side figure;imshow(subimgFlu2,[])
                
                coreFluo=mean(subimgFlu2(find(subimg3)));
                coreFluostd=std(subimgFlu2(find(subimg3)));
                perifluo=mean(subimgFlu2(find(subimg4)));
                perifluostd=std(subimgFlu2(find(subimg4)));
                if ((perifluo+0.01*perifluostd)>1.2*coreFluo)
                    continue;
                end
                
                IntNEx(k)=meanInt-meanIntBack;
                IntNExMax(k)=meanIntLocMax-meanIntBack;
                IntNExBack(k)=meanIntBack;
                IntNExBackInside(k)=meanIntBackInside;
                IntNExBackOutside(k)=meanIntBackOutside;
                intIntegratedx(k)=IntegratedNPCsignal;
                thinNE(regPropMask(k).BoundingBox(2)-8:regPropMask(k).BoundingBox(2)-8+regPropMask(k).BoundingBox(4)+16,...
                    regPropMask(k).BoundingBox(1)-8:regPropMask(k).BoundingBox(1)-8+regPropMask(k).BoundingBox(3)+16)=...
                    thinM;
                % Measure shift between Nup ring and HDEL ring precisely. Draw rays and interpolate along them
                % Find centre of Nuclear envelope
                thinMM=imdilate(thinM,strel('disk',3));
                CentreRays=[regPropMask(k).Centroid(1),regPropMask(k).Centroid(2)];
                postArray=zeros(1,90);
                for degRot=1:4:360
                    try
                        Xpts=(d/2)+cosd(degRot)*(1:3000)/100;
                        Ypts=(c/2)+sind(degRot)*(1:3000)/100;
                        VbackFlu1 = interp2(subimg,Xpts,Ypts,'spline');
                        VbackFlu2 = interp2(subimgFlu2,Xpts,Ypts,'spline');
                        Vbackg = interp2(double(thinMM),Xpts,Ypts,'bilinear');%figure;imshow(Vbackg,[])
                        VbackFlu1=VbackFlu1.*Vbackg;
                        VbackFlu2=VbackFlu2.*Vbackg;
                        
                        [maxx, post]=max(VbackFlu1);
                        [maxx2, post2]=max(VbackFlu2);
                        VbackFlu2=VbackFlu2./max(VbackFlu2(:));
                        if abs((post-post2))>200 %any shift over 2 pixels is not physical
                            continue
                        end
                        postArray(round((degRot/4))+1)=post-post2;
                    catch
                        continue
                    end
                end
                postArrayArray=[postArrayArray,mean(postArray)];
            catch
                continue % hack to avoid dealing with image boarders
            end
        end
        
        
        % measure average intensity within NE mask, restricted to in focus area
        % for entire images.
        
        IntNE(numIm)=median(IntNEx(IntNEx>0));
        IntNEMax(numIm)=median(IntNExMax(IntNExMax>0));
        IntNEBack(numIm)=median(IntNExBack(IntNExBack>0));
        IntNEBackInside(numIm)=median(IntNExBackInside(IntNExBackInside>0));
        IntNEBackOutside(numIm)=median(IntNExBackOutside(IntNExBackOutside>0));
        intIntegrated(numIm)=median(intIntegratedx(intIntegratedx>0));
        
        IntNEx=0;  IntNExBack=0;IntNExMax=0;intIntegratedx=0;
        % write NE masks for documentation purpose;
        if ~iscell(selectedImage)
            imwrite(uint16(bwmorph(thinNE,'dilate')),[selectedImage(numIm).name(1:end-4),'Mask2.tif'],'tif');
        else
            imwrite(uint16(bwmorph(thinNE,'dilate')),[pathName,selectedImage{numIm}(1:end-8),'Mask2.tif'],'tif');
        end
        postArrayArrayArray=[postArrayArrayArray,median(postArrayArray)];
    catch
        continue %control images etc need not interupt the quantification
    end
end


%create array of image file names for copy in excel files
llo=zeros(size(postArrayArrayArray,2),100)
llo=char(llo)
for oo=1:size(postArrayArrayArray,2)
    llo(oo,1:size(selectedImage(oo).name,2))=selectedImage(oo).name;
end
fileList=llo;


IntNE=IntNE';
IntNEMax=IntNEMax';
IntNEBack=IntNEBack';
IntNEBackInside=IntNEBackInside';
IntNEBackOutside=IntNEBackOutside';
postArrayArrayArray=postArrayArrayArray';
intIntegrated=intIntegrated';
fileList

B={};

for z=1
    B{z,1}='filenames';
    B{z,2}='average Intensity';
    B{z,3}='peak intensities'
    B{z,4}='background inside';
    B{z,5}='background outside'
    B{z,6}='nothing'
    B{z,7}='average integrated intensity along contours'
end

for z=1:size(fileList,1)
    B{z+1,1}=fileList(z,:);
    B{z+1,2}=IntNE(z);
    B{z+1,3}=IntNEMax(z);
    B{z+1,4}=IntNEBackInside(z);
    B{z+1,5}=IntNEBackOutside(z);
    B{z+1,6}=postArrayArrayArray(z);
    B{z+1,7}=intIntegrated(z);
end


filename = 'Intensity_NE_Output.xlsx';
xlswrite(filename,B)
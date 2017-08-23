function LocMax = LocMaximas1D( InIm, spacing )
%LocMaximas
%   identifies local maxima in images, such that they are spaced further
%   than spacing
InIm=InIm+0.000012*rand(size(InIm));

c=maxfilt2(InIm,double(spacing));
c=c(:,1);
newimg1=zeros(size(c));
try
imSameMax=c-InIm;%figure;imshow(imSameMax,[])
catch
    imSameMax=c(:,1)-InIm;%figure;imshow([InIm,300*(imSameMax==0&InIm>0)],[])

end
LocMax=find(imSameMax==0&InIm>0);

%[LocMaxImgX,LocMaxImgY]=ind2sub(size(InIm),f);

end


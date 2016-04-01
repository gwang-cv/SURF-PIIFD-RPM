function SURF_PIIFD_RPM(f1,f2)
% ========================================================================
% Robust Point Matching For Multimodal Retinal Image Registration
% Copyright(c) 2015 Gang Wang
% All Rights Reserved.
%-----------------------------------------------------------------------
% Author: Gang Wang
% Contact: gwang.cv@gmail.com
% Version: 2014-11-30 
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereQ
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
% For more information, please refer to
% Gang Wang et al., "Robust Point Matching Method for Multimodal Retinal
% Image Registration", Biomedical Signal Processing and Control, 2015, 
% Vol. 19, pp. 68-76.
%----------------------------------------------------------------------
%
% Input : (1) f1: the first image
%         (2) f2: the second image
%
% Output: fusion image
% ========================================================================

addpath('./FGT_mex');
addpath('./Utilities');
addpath('./OpenSURF_version1c');

if nargin==0
f1 = 'f1.jpg';
f2 = 'f2.jpg';
end

tic;
disp('please wait...');

if exist(f1) && exist(f2)
    im1=imread(f1);
    im2=imread(f2);
    
    im1=imresize(im1,0.3);
    im2=imresize(im2,0.4);
    
    
    % SURF
    p1=[]; p2=[];
    Options.upright=true;
    Options.tresh=0.00001;
    Ipts1=OpenSurf(im1,Options);
    Ipts2=OpenSurf(im2,Options);
    for iss=1:length(Ipts1)
        p1s=[Ipts1((iss)).x ,Ipts1((iss)).y];
        p1=[p1;p1s];
    end
    for iss=1:length(Ipts2)
        p2s=[Ipts2((iss)).x ,Ipts2((iss)).y];
        p2=[p2;p2s];
    end
    p1=round(p1);
    p2=round(p2);
    %proprocessing
    if size(im1,3)>1
        im1=im1(:,:,2);
    end
    if size(im2,3)>1
        im2=im2(:,:,2);
    end
    im1=double(im1);
    im2=double(im2);
    im1=im1/max(im1(:));
    im2=im2/max(im2(:));
    I1=im1;
    I2=im2;
    I1 = I1*255;
    I2 = I2*255;
    
    %----- remove boundry's points ------%
    [msk1,msk2] = rr_msk(I1,I2,15);
    p1= rr_removeboundarypoint(p1,msk1,2);
    p2= rr_removeboundarypoint(p2,msk2,2);
    cols1=p1(:,1);
    cols2=p2(:,1);
    rws1=p1(:,2);
    rws2=p2(:,2);
    % show corner points
    figure;
    imshow(im1); hold on, plot(p1(:,1),p1(:,2),'r*');
    figure;
    imshow(im2); hold on, plot(p2(:,1),p2(:,2),'ro');
    
    match=rr_desmatch(I1,I2,cols1,cols2,rws1,rws2);   
    loc1 = match(:,1:4);
    loc2 = match(:,5:8);
    showmatch(im1,im2,loc1,loc2,0);
    
    %----- Roubst Point Matching with outliers ------%
    x1=match(:,1:2);
    x2=match(:,5:6);
    [ys,indx]=RPM(x1,x2);
    loc1=match(indx,1:4);
    loc2=match(indx,5:8);
    loc1(:,1:2) = cpcorr(loc1(:,1:2),loc2(:,1:2),im1,im2);
    showmatch(im1,im2,loc1,loc2,0);
    
    %--------Estimate Transformation ----------%
    pointNum=size(loc1,1);
    if pointNum >=3
        % Affine
        t_fundus = cp2tform(loc1(:,1:2),loc2(:,1:2),'affine');
        I1_c = imtransform(im1,t_fundus,'XData',[1 size(im2,2)], 'YData',[1 size(im2,1)]);
        [I1_c,I2_c] = rr_imagesize(I1_c,im2);
        a1=figure;imshow(I1_c+I2_c,[]); title(['fusion image (affine)']);% intensity
        disp('Affine transformation is done.');
        if pointNum >=6
            % 2nd Polynomial
            t_fundus = cp2tform(loc1(:,1:2),loc2(:,1:2),'polynomial',2);
            I1_c = imtransform(im1,t_fundus,'XData',[1 size(im2,2)], 'YData',[1 size(im2,1)]);
            [I1_c,I2_c] = rr_imagesize(I1_c,im2);
            po1=figure;imshow(I1_c+I2_c,[]);title(['fusion image (Polynomial)']);
            disp('2nd Polynomial transformation is done.');
        end
    else
        fprintf('need at least 3 points !')
    end
else
    fprintf('The files can not be found !')
end
toc;

%% utility functions
function showmatch(I1_p,I2_p,loc1,loc2,s,fileNum,label)

I1 = I1_p(s+1:end-s,s+1:end-s);
I2 = I2_p(s+1:end-s,s+1:end-s);
loc1 = loc1-s;
loc2 = loc2-s;

im3 = rr_appendimages(I1,I2);

figure,imshow(im3,[])
title(['matched points (' num2str(size(loc1,1)) ')']);
hold on
cols = size(I1,2);

for i=1:size(loc1,1)
    line([loc1(i,1) loc2(i,1)+cols],[loc1(i,2) loc2(i,2)], 'Color', 'y');
    plot(loc1(i,1),loc1(i,2),'g.')
    plot(loc1(i,1),loc1(i,2),'go')
    plot(loc2(i,1)+cols,loc2(i,2),'g.')
    plot(loc2(i,1)+cols,loc2(i,2),'go')
end













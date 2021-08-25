function [pkcentres] = findpeak(zim,param)
%%%%%%%%%%%%%%%%%%
% Name: findpeak_timeTotrack
% Purpoise: find the peaks positions inside the cell and prepare the file to give to the
% track function
%
% INPUT:
% zim: images organised in a cell file (see read_stacks)
% param: parameters for find the peaks
%    param.thrfpeak=12.5;  %% multiplicative factor for trh find peak function
%    param.pnoise=1;       %% noise in peakfind function
%    param.psize=4;        %% size group of pixel in peakfind function
%    param.pgauss=7;       %% size for gaussian fit peak size in centfind
% OUTPUT: CHANGED
% file(table) organize like
%                    (x)      (y)      (t)
% ;     pos = 3.60000      5.00000      0.00000
% ;           15.1000      22.6000      0.00000
% ;           4.10000      5.50000      1.00000
% ;           15.9000      20.7000      2.00000
% ;           6.20000      4.30000      2.00000
% to give to the track function
%
% function developed by Alessia Lepore- El Karoui lab 2017
% other functions used here(see documentation):
% pkfind
% bpass
% centfind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Filter image and find peaks
imfilt=zeros(size(zim));
pkpos=[];
for i = 1:param.nstacks
    imfilt(:,:,i) = bpass(zim(:,:,i),param.pnoise,param.psize); % Image filter (bandpass)
    thr = mean2(imfilt(:,:,i));
    peaks = pkfnd(imfilt(:,:,i),thr*param.thrfpeak,param.psize); % Find intensity maxima
    pkpos = [pkpos; peaks,repmat(i,size(peaks,1),1)];
end

% Fit gaussian for sub-pixel localisation
pkcentres = zeros(size(pkpos,1),7);
for i = 1:max(pkpos(:,3))
    clc
    disp(['Peak finder processing frame ' num2str(i) '/' num2str(max(pkpos(:,3)))]);
    [~, pkcentres(pkpos(:,3)==i,:)] = centfind(imfilt(:,:,i),pkpos(pkpos(:,3)==i,1:2),param.pgauss,0,['Gaussian','interactive']);
end

pkcentres(:,8) = pkpos(:,3); % Add frame numbers to output



%% Show localisations on raw and filtered images

% framesToShow = 1:param.disp_rate:param.nstacks;
% c = 0;
% figure('Color','white')
% for i = framesToShow
%     c = c+1;
%     if ~sum(pkpos(:,3)==i)
%         subplot(ceil(length(framesToShow)/2),2,c)
%         imagesc(zim(:,:,i))
%         axis equal
%         colormap('gray')
%         title(['Frame ' num2str(i)]);
%         hold on
%         plot(pkcentres(pkcentres(:,8)==i,5),pkcentres(pkcentres(:,8)==i,6),'.r','MarkerSize',6)
%     else
%         j = i;
%         while ~sum(pkpos(:,3)==j)
%             j = j + 1;
%         end
%         subplot(ceil(length(framesToShow)/2),2,c)
%         imagesc(zim(:,:,j))
%         axis equal
%         colormap('gray')
%         title(['Frame ' num2str(j)]);
%         hold on
%         plot(pkcentres(pkcentres(:,8)==j,5),pkcentres(pkcentres(:,8)==j,6),'.r','MarkerSize',6)
%     end
% end
% hold off
% 
% 
% 
% framesToShow = 1:param.disp_rate:param.nstacks;
% c = 0;
% figure('Color','white')
% for i = framesToShow
%     c = c+1;
%     if ~sum(pkpos(:,3)==i)
%         subplot(ceil(length(framesToShow)/2),2,c)
%         imagesc(imfilt(:,:,i))
%         axis equal
%         colormap('gray')
%         title(['Frame ' num2str(i)]);
%         hold on
%         plot(pkcentres(pkcentres(:,8)==i,5),pkcentres(pkcentres(:,8)==i,6),'.r','MarkerSize',6)
%     else
%         j = i;
%         while ~sum(pkpos(:,3)==j)
%             j = j + 1;
%         end
%         subplot(ceil(length(framesToShow)/2),2,c)
%         imagesc(imfilt(:,:,j))
%         axis equal
%         colormap('gray')
%         title(['Frame ' num2str(j)]);
%         hold on
%         plot(pkcentres(pkcentres(:,8)==j,5),pkcentres(pkcentres(:,8)==j,6),'.r','MarkerSize',6)
%     end
% end
% hold off




end
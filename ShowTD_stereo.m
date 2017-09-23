function  ShowTD_stereo(varargin)
%   Shows a video of Temporal Difference (TD) stereo events
%   function  ShowTD_stereo_simple(TD,FPS,output_filename)
% TAKES IN:
%   'TD'
%       A struct of TD events with format:
%           TD(:,1) = event timestamps in microseconds
%           TD(:,2) = pixel X locations
%           TD(:,3) = pixel Y locations
%           TD(:,4) = event polarity
%           TD(:,5) = left or right label left is 0 right is 1
%           TD(:,6) = disparity
%
%   'FPS'
%       Defaults to 24FPS, which is a Time Per
%       Frame (TPF) of 1/24 seconds.
%
%    'Filename'
%       Record the video if filename is not empty
% written by Xie zhen -May 2017
% zjutxz@hotmail.com

timeconst = 1e-6;
TD = varargin{1};

%FPS is 1/TPF
if nargin > 1
    if isempty(varargin{2})
        TPF = 1/(50*timeconst);
    else
        TPF = 1/(varargin{2}*timeconst);
    end
else
    TPF = 1/(50*timeconst);
end

if nargin > 2
    if isempty(varargin{3})
        aviname = 'stereo_input.avi';
    else
        aviname = strcat(varargin{3},'.avi');
        
    end
    % to create the video
    aviobj=VideoWriter(aviname,'Motion JPEG AVI');
    aviobj.Quality = 100;
    open(aviobj);
    % set the figure
    f = figure(1);
    axis off;set(f,'menubar','none','toolbar','none');
    hold on
else
    disp('Do not record the video');
end


%***** modified when use different sensor*******
% 180*240 for DAVIS
weight = 240;
height = 180;

ImageBack_LR = zeros(height,weight,3);

% make sure p > 0
if min(TD(:,4))==0
    TD(:,4) = TD(:,4)+1;
end

% define color map
cc = hsv(double(max(TD(:,4))));
Image_L = ImageBack_LR;
Image_R = ImageBack_LR;

estimate_interval_number = 0;

for  i = 3: length(TD)
    if (estimate_interval_number ~=  floor((TD(i,1)-TD(1,1))/TPF) )
        estimate_interval_number = floor((TD(i,1)-TD(1,1))/TPF);
        Image = [Image_L zeros(size(Image_R,1),2,3)+255 Image_R];
        imshow(Image);
        drawnow();        
        if nargin > 2
            m=getframe(f);
            writeVideo(aviobj,m);
        end
        Image_L = ImageBack_LR;
        Image_R = ImageBack_LR;
    else
        if isequal(TD(i,5), 0)
            Image_L(TD(i,3), TD(i,2), :) = cc(TD(i,4),:);
        elseif isequal(TD(i,5), 1)
            Image_R(TD(i,3), TD(i,2), :) = cc(TD(i,4),:);
        end
    end
end

end

function ShowTD_depth(varargin)
% ShowTD_depth(TD, FPS, Filename)
%   Shows a video of Temporal Difference (TD) events depth and save a video
%   All arguments except TD are optional.
%
% TAKES IN:
%   'TD'
%       A struct of TD events with format:
%           TD.x =  pixel X locations
%           TD.y =  pixel Y locations
%           TD.p =  event polarity
%           TD.ts = event timestamps in microseconds
%           TD.disparity = disparity of events in pixel
%
%
%    'Filename'
%       Record the video if filename is not empty
%     
% written by Zhen Xie - Feb 2017


% close all
timeconst = 1e-6;
TD_1 = varargin{1};
%   'FPS'
%       Defaults to 48FPS, which is a Time Per
%       Frame (TPF) of 1/24 seconds.
FPS = 24;
%     'FL' overlap
%       Frame Length (FL) is an optional arguments specifying the time-span
%       of data to show per frame as a fraction of TPF. Defaults to TPF seconds. If FL<1,
%       then not all data in a sequence will be shown. If FL>1 then some
%       data will be repeated in subsequent frames.
Overlap = 1;

Tmin = 1;
Tmax = length(TD_1.ts);

if nargin > 1
    if isempty(varargin{2})
        aviname = 'result.avi';
    else
        aviname = strcat(varargin{2},'.avi');
    end
    % to create the video
    aviobj=VideoWriter(aviname,'Motion JPEG AVI');
    aviobj.Quality = 100;
    open(aviobj);
    % set the figure
    f = figure(1);
    axis off;%set(f,'menubar','none','toolbar','none');
    hold on
else
    disp('Do not record the video');
end

FrameLength = 1/(FPS*timeconst);
t1_1 = TD_1.ts(Tmin) + FrameLength;
t2_1 = TD_1.ts(Tmin) + FrameLength*Overlap;

% for DAVIS
ImageBack = zeros(180,240,3);
if (min(TD_1.x) == 0 || min(TD_1.y) == 0)
    TD_1.x = TD_1.x + 1;
    TD_1.y = TD_1.y + 1;
end

i = Tmin;
%the range of the depth
min_d = 0.6;   
max_d = 5;   
cc = jet(double(length(min_d:0.2:max_d))); 
labels = min_d:0.2:max_d;
for i = length(labels):-1:1
    string_num_x{i} = num2str(labels(i));
end

Image1= ImageBack;
Image2 = ImageBack;
k=1;

while (i<Tmax)
    j=i;
    while ((j<Tmax) &(TD_1.ts(j) < t2_1))
        if (TD_1.p(j) ~= inf) & ~isnan(TD_1.p(j)) & ~isnan(TD_1.disparity_gt(j)) & TD_1.disparity_gt(j) ~= 0   
             if TD_1.p(j)<=6
                Image1(TD_1.y(j), TD_1.x(j), :) = cc(1,:);
             else
                cc_index = length(min_d:0.2:max_d)-round((250*0.12/TD_1.p(j)-min_d)/0.2);
                Image1(TD_1.y(j), TD_1.x(j), :) = cc(cc_index,:);
             end
             
             if TD_1.disparity_gt(j)>=5   % 252*0.12/
                Image2(TD_1.y(j), TD_1.x(j), :) = cc(1,:);
             elseif TD_1.disparity_gt(j) > min_d
                cc_index =length(min_d:0.2:max_d)-round((252*0.12/TD_1.disparity_gt(j)-min_d)/0.2);
                Image2(TD_1.y(j), TD_1.x(j), :) = cc(cc_index,:);
             else
                Image2(TD_1.y(j), TD_1.x(j), :) = cc(length(min_d:0.2:max_d),:);
             end
        end
        j = j+1;
    end
    while ((TD_1.ts(i) < t1_1) & (i<Tmax))
        i = i+1;
    end
    
    f = figure(2);
     axis off;
     colormap(flipud(cc));
     lcolorbar(string_num_x);
     Image = Image1;
     imshow(Image);
     drawnow();
     
     if nargin > 2
         m=getframe(f);
         writeVideo(aviobj,m);
     end
     
 
    t2_1 = t1_1 + FrameLength*Overlap;
    t1_1 = t1_1 + FrameLength;
    
    Image1 = ImageBack;
    Image2 = ImageBack;
    k=k+1;
end
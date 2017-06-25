# EMP
event-based stereo matching algorithm
This README would normally document whatever steps are necessary to get your application up and running.
What is this repository for?
The event-based stereo uses event based camera such as DAVIS ATIS.
The input is each event e(ts,x,y,p,l,Gt_d), ts is timestamps, x and y are the image coodrinates of the events,p is the polarity, l is the left/right label(0 is left and 1 is right),Gt_d is the ground truth of disparity.
The output is the disparity value for each event. 
This work is showed detailed in Xie Zhen, Garrick Orchard,"Asynchronous event based stereo matching by message passing algorithm" 2017 If you use this work in an academic context, please cite the paper.

Any questions
Please Contact Harry xie xzzjut@163.com

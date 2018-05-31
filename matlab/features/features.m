%features  class to contain and keep tracked frames trough a video.
%
% Syntax:
%    obj = features(p_frames,N_features,p_camera_model)
%
% In:
%    p_frames - a cell array containing the path to the frames
%    N_features - the maximum number of features to track.
%    p_camera - model is the camera model (optional) 
%
% Description:
%   Track salient features in a video using the good features to track
%   algorithm. store the tracks and return them.
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

classdef features < handle
    %Class for feature handling.
    
    properties
        frames; %Paths of frame images.
        im_size;%Size of the image.
        camera_model; %Parameters of camera.
        t; %Current frame.
        N; %Length of feature vector.
        coord; % Matrix containing the coordinates((2*N_frames) x N ).
        event; %(sparse)Matrix containing death an birth of feature tracks 0-alive, 1-dead, 2-birth.
    end
    
    methods
        
        function obj = features(p_frames,N_features,p_camera_model)
            %Class initializer, where p_frames is a cell array containing
            %the path to the frames, N_features are the maximum number of
            %features to track, and p_camera model is the camera model
            %(optional).
            
            obj.frames=p_frames;
            obj.camera_model=p_camera_model;
            obj.N=N_features;
            obj.coord=[];
            obj.t=1;
            obj.event=sparse([]);
        end
        
        function S = draw_all(obj)
            %Draw all features with a tail of 10 frames, usefull for
            %spotting jumps among the tracks.
            
            for i=1:length(obj.frames)-1
                obj.draw_tracks(i,min(i,10),0);
                drawnow();
            end
        end
        
        function S= refine_track(obj,im)
            %Refine the track of the current frame by using the 7 point
            %algorithm as described in Hartley-Zissermans book.
            
            %Read the coordinates tracked for the current and previous
            %frame.
            [S,status,~] = obj.read_tracks(obj.t,1,0);
            %Find features inside the frame.
            in_frame=(~(S(:,3)<=0 |  S(:,3)>obj.im_size(1))) & (~(S(:,4)<=0 |  S(:,4)>obj.im_size(2)));
            P=S(status(:,1) & in_frame,:);
            %Find inliers within the current frame.
            inl=ransac_7_point(P,1,0.75);
            
            %Track inlier indexes back to original array.
            label=1:obj.N;
            inliers=label(status(:,1) & in_frame);
            outliers=(union(label(~status(:,1)),setdiff(inliers,inliers(inl))));
            
            
            %Measure euclidean difference between frames.
            dife=(sum((S(:,3:4)-S(:,1:2)).^2,2).^0.5);
            
            %stop big gaps
            outliers=union(outliers,label(dife>40));
            
            %Brand outliers as such.
            obj.event(outliers,obj.t-1)=1;
            obj.event(outliers,obj.t)=2;
            
           
            
           
            
            %Diagnosis with maximum distance.
            final_inliers=setdiff(label,outliers);
            diag=max(dife(final_inliers))
            1;
            
             %Restart tracks.
            obj.reborn(im);
            
        end
        
        function S = draw_tracks(obj,t,Max_back,Max_forward)
            %Draw all tracks alive at a frame t with a backward tail of
            %Max_back and forward Max_forward.
            
            %Read coordinates for the given frames.
            [tracks,mask,tstamp,eve]=read_tracks(obj,t,Max_back,Max_forward);
            
            %draw frame to initialize figure.
            obj.draw(tstamp(1))
            hold on
            
            %For every frame.
            for i =1:length(tstamp)
                %Draw current frame.
                obj.draw(tstamp(i))
                %For every feature.
                for j=1:obj.N
                    
                    %Find alive features.
                    vecx=(1:2:2*i-1);
                    vecy=(2:2:2*i);
                    x=tracks(j, vecx( mask(j,1:i)==1));
                    y=tracks(j, vecy( mask(j,1:i)==1));
                    %draw line
                    line(x,y,'Color','g')
                    
                end
            end
        end
        
        
        function S=track_all(obj)
            %track all frames from first one till the end, this includes
            %refining tracks based on stereo triangulation.
            
            %track first frame.
            im=obj.track_first();
            for i=1:length(obj.frames)-1
                %track next frame.
                im=obj.track_next(im);
                %refine current track.
                %obj.refine_track(im);
                %draw features
                obj.draw(obj.t);
                
                %Draw tracks (slow)
                if obj.t==140 && 0
                    obj.draw_tracks(obj.t,obj.t-1,0);
                    keyboard;
                end
                
                drawnow();
            end
            %Delete last event since its false initialized in the track
            %next method.
            obj.event(:,end)=[];
        end
        
        function [S,status,tst,eve] = read_tracks(obj,t,Max_back,Max_forward)
            %Return all active tracks at time t for the range indicated in
            %the parameters.
            
            %Find indexes for the frames.
            tst=max(1,t-Max_back):min(length(obj.frames),(t+Max_forward));
            
            %Index of central frame.
            t_ind=t-max(1,t-Max_back);
            
            %Index relative to current selection.
            ind=tst-t;
            
            %Read track status.
            eve=obj.event(:,tst);
            status=ones(size(eve));
            
            %From the middle propagate alive status back.
            for i=t_ind:-1:1
                status(:,i)=status(:,i+1)&eve(:,i)~=1;
            end
            
            %From the middle propagate alive status forward
            for i=t_ind+2:1:length(tst)
                status(:,i)=status(:,i-1)&eve(:,i-1)~=1;
            end
            
            %Read coordinates of alive tracks.
            S=obj.coord(:,(tst(1))*2-1:(tst(end))*2);
            S(:,1:2:end)=S(:,1:2:end).* status;
            S(:,2:2:end)=S(:,2:2:end).* status;
            status=status==1;
        end
        
        function S=draw(obj,t)
            %Draw a frame and the features tracked on it.
            im = rgb2gray(imread(obj.frames{t}));
            colormap('gray')
            if ~isempty(obj.camera_model)
                im=undistortImage(im);
            end
            imagesc(im);
            hold on
            co=obj.coord(:,(2*t)-1:2*t);
            ind=obj.event(:,t);
            plot(co(ind==0,1),co(ind==0,2),'ro');
            plot(co(ind==2,1),co(ind==2,2),'go');
            hold off
            
        end
        
        
        function im=track_first(obj)
            %Track first frame.
            
            %mark all as born (initializing size)
            obj.event(1:obj.N,1)=2;
            obj.coord(1:obj.N,1:2)=0;
            %track first frame
            im = rgb2gray(imread(obj.frames{obj.t}));
            if ~isempty(obj.camera_model)
                im=undistortImage(im);
            end
            mask=ones(size(im));
            obj.im_size=size(im);
            %find features
            pts_add = cv.goodFeaturesToTrack(im,'MaxCorners',obj.N,'Mask',mask,'MinDistance',20);
            
            %fill status and data with results of feature search
            index_born=1:length(pts_add);
            index_postponed=(length(pts_add)+1):obj.N;
            obj.event(index_postponed,obj.t+1)=2;
            obj.event(index_postponed,obj.t)=1;
            obj.coord(index_born,(obj.t)*2-1:(obj.t)*2)=cell2mat(pts_add');
        end
        
        
        function im = track_next(obj,im_prev)
            %track next frame with the current features as initial points.
            
            %move to next frame
            obj.t=obj.t+1;
            
            %read Image (undistort if a camera object is provided)(to be changed to better calibration format)
            im = rgb2gray(imread(obj.frames{obj.t}));
            if ~isempty(obj.camera_model)
                im=undistortImage(im);
            end
            
            %find alive points to continue track
            index_alive=obj.event(:,obj.t)~=(1);
            pts=obj.coord(index_alive,(obj.t-1)*2-1:(obj.t-1)*2);
            [pts_new,stat] = cv.calcOpticalFlowPyrLK(im_prev, im, pts,'WinSize',[15 15]);
            
            %kill tracks that are interrumpted and mark new ones as born
            obj.event((stat==0),obj.t-1)=1;
            obj.event((stat==0),obj.t)=2;
            %index_alive=index_alive(status==1);
            obj.coord(index_alive,obj.t*2-1:obj.t*2)=cell2mat(pts_new');
            index_born=obj.event(:,obj.t)==(2);
            
            %create mask
            mask=ones(size(im));
            for i=1:size(pts,1)
                imask = 5;
                indj = pts(i,1)+(-imask:imask);
                indj = min(indj,size(im,2)); indj = max(indj,1);
                indi = pts(i,2)+(-imask:imask);
                indi = min(indi,size(im,1)); indi = max(indi,1);
                mask(round(indi),round(indj)) = 0;
            end
            obj.event(:,obj.t+1)=0;
            if sum(index_born)>0
                %track points to be born
                pts_add = cv.goodFeaturesToTrack(im,'MaxCorners',sum(index_born),'Mask',mask,'MinDistance',5);
                
                
                cumu=cumsum(index_born);
                %fill status and data with results of feature search
                index_newborn= index_born & (cumu<=length(pts_add));
                index_postponed= index_born & (cumu>length(pts_add));
                obj.event(index_postponed,obj.t+1)=2;
                obj.event(index_postponed,obj.t)=1;
                obj.event(index_newborn,obj.t+1)=0;
                obj.coord((index_newborn),obj.t*2-1:obj.t*2)=cell2mat(pts_add');
            end
            
            
        end
        
        function S = reborn(obj,im)
            %restart features branded to be born on current frame
            
            index_born=obj.event(:,obj.t)==(2);
            index_alive=obj.event(:,obj.t)==(0);
            pts=obj.coord(index_alive,(obj.t)*2-1:(obj.t)*2);
            
            %create mask
            mask=ones(obj.im_size);
            
            %For every frame
            for i=1:size(pts,1)
                %Size of square mask.
                imask = 5;
                
                %Find indexes corresponding to mask
                indj = pts(i,1)+(-imask:imask);
                indj = min(indj,size(mask,2)); indj = max(indj,1);
                indi = pts(i,2)+(-imask:imask);
                indi = min(indi,size(mask,1)); indi = max(indi,1);
                
                %Create mask.
                mask(round(indi),round(indj)) = 0;
            end
            
            
            obj.event(:,obj.t+1)=0;
            if sum(index_born)>0
                %track points to be born
                pts_add = cv.goodFeaturesToTrack(im,'MaxCorners',sum(index_born),'Mask',mask,'MinDistance',10);
                
                
                cumu=cumsum(index_born);
                %fill status and data with results of feature search
                index_newborn= index_born & (cumu<=length(pts_add));
                index_postponed= index_born & (cumu>length(pts_add));
                obj.event(index_postponed,obj.t+1)=2;
                obj.event(index_postponed,obj.t)=1;
                obj.event(index_newborn,obj.t+1)=0;
                obj.event(index_newborn,obj.t)=2;
                obj.coord((index_newborn),obj.t*2-1:obj.t*2)=cell2mat(pts_add');
            end
        end
        
    end
    
end
%Example of selfcalibration using gyro and image data.
%
% Copyright (C) Santiago Cortés
%
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

%Enable frame by frame visualization of visual measurements.
visualize=false;

%Seed random number generator
rng(0)

%define indeces for state and covariance matrices
%c=fx,fy,cx,cy.
%p is position in R3.
%v is velocity in R3
%q is orientation as a quaternion.
%w is rotation rate in R3.
%w is the gyroscope data, in R3.
%z is the feature vector positions (N by 3).
%x= c p v q w z;

inde.c=1:6;
inde.p=inde.c(end)+1:inde.c(end)+3;
inde.v=inde.p(end)+1:inde.p(end)+3;
inde.pv=[inde.p inde.v];
inde.q=inde.v(end)+1:inde.v(end)+4;
inde.w=inde.q(end)+1:inde.q(end)+3;

%number of features
N_feat=6*9;

%index for 3d locations
inde.z=inde.w(end)+1:inde.w(end)+N_feat*3;
init=zeros(3,N_feat);

%initial scale parameter
siz=1;
%initial position
init(3,:)=siz*(randn([N_feat,1])*1+10);
init(1,:)=siz*randn([N_feat,1])*1;
init(2,:)=siz*randn([N_feat,1])*1;

%% Load data
folder='../data/';
datapath_calib = 'cards2';

%Load data dump if it exists, if not, create and save data dump.
if exist([folder 'data_' datapath_calib '_' int2str(N_feat) '_features'  '.mat'],'file')
    load([folder 'data_' datapath_calib '_' int2str(N_feat) '_features'  '.mat'])
    %If no datadump exists perform the tracking (requires mexopencv)
else
    datapath_patt =[folder 'calib'] ;
    [cameraParams,cameraMatrix] = calibrateCameraMatlab(datapath_patt,1);
    [data,frames] = read_camera_data([folder datapath_calib],true);
    feat=features(frames,N_feat,[]);
    feat.track_all();
    wb= [0 0 0];
    save([folder 'data_' datapath_calib '_' int2str(N_feat) '_features'  '.mat'],'feat','cameraMatrix','wb','data','frames','N_feat')
end
if visualize
    [data,frames] = read_camera_data([folder datapath_calib],false);
end

%% align timestamps.
ind = (data(:,2)==22);
data(ind,1) = data(ind,4)-0.1;

% Sort.
data = sortrows(data);

% Cut away non-valid parts.
data = data(max(find(data(:,2)==22,1),find(data(:,2)==4,1)): ...
    min(find(data(:,2)==22,1,'last'),find(data(:,2)==4,1,'last')),:);

% Start at zero.
data(:,1) = data(:,1) - min(data(:,1));
eve=[];

%%
%initialize state vector.
m(inde.q)=[1 0 0 0];
m(inde.w)=[0 0 0];
m(inde.p)=[0 0 0];
m(inde.v)=[0 0 0];
m(inde.c)=[800 800 240 340 0 0];
m(inde.z)=init(:);
%inititialize time.
gyroT=0;
cameraT=0;
ind=0;
%% Dynamic model
%wiener velocity model, random walk in acceleration (smooth motion prior).

%noise parameter for wiener velocity model (acceleration units), is kept
%relative to the scene scale .
qvel=2;
w_var=[1 1 1]*0.01^2;
%square root of diagonal of Q matrix (process noise matrix), in the same units as
%the corresponding state variable.
delT=0.1;
di=m*0;
di(inde.w)=0;
di(inde.c)=0;
di(inde.z)=0;
di(inde.q)=0;

%%square root of diagonal of initial P matrix (covariance matrix), in the same units as
%the corresponding state variable.
pi=m;
pi(inde.p)=0;
pi(inde.v)=0;
pi(inde.q)=0;
pi(inde.w)=0;
pi(inde.c)=[8 8 8 8 0.1 0.1 ]*10;
pi(inde.z)=1;
%Covariance and "square root" of covariance matrix.
P=diag(pi.^2);
Skk=diag(pi);

%create base for Q matrix
Qs=diag(di.^2);
Qs(inde.pv,inde.pv)=wienVel(qvel,delT);

%measurement noise matrix (in pixels squared).
Rcam=eye(N_feat*2)*12^2;

%loop control variables
exit=false;
ki=0;
k=0;

for ind=1:size(data,1)
    
    %read gyro data
    if  data(ind,2)==4
        %advance control step.
        k=k+1;
        
        %read gyro data and remove measured bias.
        w=data(ind,3:5)-wb;
        
        %fix handedness and orientation of gyro measuremet.
        w=[w(1) -w(2) -w(3)];
        
        %time step
        delT=data(ind,1)-gyroT;
        
        %scale process noise according to timestep
        Q=Qs;
        Q(inde.pv,inde.pv)=wienVel(qvel,delT);
        Q(inde.q,inde.q)=quat_cov(m(inde.q),w_var,delT);
        gyroT=data(ind,1);
        
        %matrix form of quaternion operator
        crow=[0 w(3) -w(2);-w(3) 0 w(1);w(2) -w(1) 0];
        ome=[0 -w;w' crow ];
        ex_ome=expm(delT*ome/2);
        
        %create linear process matrix
        A=eye(length(m));
        A(inde.pv,inde.pv)=[eye(3) eye(3)*delT;zeros(3) eye(3)];
        A(inde.q,inde.q)=ex_ome;
        A(inde.w,inde.w)=diag(w);
        
        %make kalman prediction step
        [m,P] = kf_predict(m',P,A,Q);
        
        %keep quaternion normalized
        m=m';
        m(inde.q)= m(inde.q)/norm( m(inde.q));
             
        %store result.
        qte(k,:)=m(inde.q);
        MEAN(k,:)=m;       
        VAR(k,:)=diag(P);
        
        
        if visualize
            %project points for visualization.
            points=m(inde.z);
            points=reshape(points,3,N_feat);
            points=points(:,eve~=1);
            points=[points;ones(1,size(points,2))];
            
            %draw projected points
            c=m(inde.c);
            K=[c(1) 0 c(3);0 c(2) c(4); 0 0 1];            
            R=quat2rmat(m(inde.q));
            p=m(inde.p);
            E=[R' -R'*p'];
            Y=K*E*points;
            S=Y(1:2,:)./repmat(Y(3,:),2,1);
            figure(1)
            hold on
            scatter(S(1,:),S(2,:),'cx');
            drawnow;
        end
        
    end
    
    %Image data
    if data(ind,2)==22
        ki=ki+1;
        if visualize
            %Draw image in the background
            figure(1)
            im=rgb2gray(imread(frames{ki}));         
            hold off
            imagesc(repmat(im,[1 1 3]));
            hold on
        end

        %read feature position
        [coord,status,tst,eve]=read_tracks(feat,ki,0,0);
        resind1=full(eve==1 | eve==2);
        resind2=full(eve==1);
        inl_state=(ones(size(m))==1)';
        inl_meas=ones(numel(inde.z)/3*2,1)==1;
        
        %reset 3d state of lost features
        if sum(resind1)>0            
            index_to_reset_3d(1:3:N_feat*3)=resind1;
            index_to_reset_3d(2:3:N_feat*3)=resind1;
            index_to_reset_3d(3:3:N_feat*3)=resind1;
            
            index_to_ignore_3d(1:3:N_feat*3)=resind2;
            index_to_ignore_3d(2:3:N_feat*3)=resind2;
            index_to_ignore_3d(3:3:N_feat*3)=resind2;
            
            clear eve1
            eve1(1,:)=eve==1;
            eve1(2,:)=eve==1;
            eve1=eve1';
            eve1=eve1(:);            
            eve2(1:2:N_feat*2)=eve==2;
            eve2(2:2:N_feat*2)=eve==2;
            
            
            %find center of point cluster
            q=m(inde.q);          
            R=quat2rmat(q);
            center=m(inde.p);
            co=sum(resind1);
            init=[];
            
            %initialize randomly around the center of the cluster
            init(3,:)=siz*(randn([co,1])*0.1);
            init(1,:)=siz*randn([co,1])*1;
            init(2,:)=siz*randn([co,1])*1;            
            vec=[];    
            points=m(inde.z);
            points=reshape(points,3,N_feat);
            points=points(:,eve~=1);           
            cloud_center=mean(points,2);
            distance_to_mean=norm(cloud_center- center');            
            for i=1:co
                vec(:,i)=cloud_center+init(:,i)*1;                
            end
            
            %clear covariance cross terms and initialize varaince of new
            %positions.
            m(inde.z(index_to_reset_3d))=vec(:);
            P(inde.z(index_to_reset_3d),:)=0;
            P(:,inde.z(index_to_reset_3d))=0;
            P(inde.z(index_to_reset_3d),inde.z(index_to_reset_3d))=(eye(sum(index_to_reset_3d))*10).^2;            
            inl_state(inde.z(index_to_ignore_3d))=0;
            inl_meas=full(eve1~=1)';
        end
              
        z=coord(:,end-1:end);
        z=z(:);
        z=z(inl_meas);
        %perform EKF update only on inlier features
        [mn,Pn] = ekf_update1(m(inl_state)',P(inl_state,inl_state),z,@(x) jac_h(x,inde),Rcam(inl_meas,inl_meas),@(x) h(x,inde));       
        m(inl_state)=(mn)';
        P(inl_state,inl_state)=Pn;
        %renormalize quaternion
        m(inde.q)= m(inde.q)/norm( m(inde.q));
        
        %show location of projected points
        if visualize
            points=m(inde.z);
            points=reshape(points,3,N_feat);
            points=points(:,eve~=1);
            points=[points;ones(1,size(points,2))];
            c=m(inde.c);
            Kma=[c(1) 0 c(3);0 c(2) c(4); 0 0 1];
            Ro=quat2rmat(m(inde.q));
            p=m(inde.p);
            E=[Ro' -Ro'*p'];
            Y=Kma*E*points;
            imp=Y(1:2,:)./repmat(Y(3,:),2,1);            
            scatter(imp(1,:),imp(2,:),'rx');
            drawnow;
        end
        
        % to print progress uncomment this lines
       
        [m(inde.c(1:4))' diag(P(inde.c(1:4),inde.c(1:4)))]
        %[ [m(inde.z(1:3)), 0]' [diag(P(inde.z(1:3),inde.z(1:3))); 0]  ]
        %[m(inde.c(5:6))' diag(P(inde.c(5:6),inde.c(5:6)))]
 
        
        p(ki,:)=m(inde.p);
        men(ki,:)=m;       
    end
    
    %if ki==2
    % S=h(m',m(inde.q),inde);
    % sca=S(end);
    
    %end
    
end


%% draw variable evolution

    tg=data(data(:,2)==4,1);
    camMat=[575 0 240;0 575 320;0 0 1];
    close all
    %    addpath('matlab2tikz');
    cam=men(:,inde.c);
    v=1;%ones(size(cam,1),1);
    gt=[v*cameraMatrix(1,1) v*cameraMatrix(2,2) v*cameraMatrix(1,3) v*cameraMatrix(2,3)];
    
    
    
    t=tg;
    
    figure
    % plot(cam(:,1));
    hold on
    % plot(gt(:,1),'k--');
    % axis([0 inf 400 800]);
    
    sampl=1:10:length(tg);
    s=VAR(sampl,inde.c(1)).^.5;
    x=MEAN(sampl,inde.c(1));
    t=tg(sampl);
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    plot([t(1) t(end)],[1 1]*gt(1),'--k')
    set(gca,'Layer','Top')
    axis([0 inf 400 1000]);
    
    ylabel('$f_x$ (px)')
    %     matlab2tikz('./../paper/Oleaf/fig/real_focal_x.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)
    
    %subplot(612)
    figure
    % plot(cam(:,2));
    hold on
    % plot(gt(:,2),'k--');
    % axis([0 inf 525 650]);
    
    s=VAR(sampl,inde.c(2)).^.5;
    x=MEAN(sampl,inde.c(2));
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    plot([t(1) t(end)],[1 1]*gt(2),'--k')
    set(gca,'Layer','Top')
    axis([0 inf 400 1000]);
    
    ylabel('$f_y$ (px)')
    
    %     matlab2tikz('./../paper/Oleaf/fig/real_focal_y.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)
    
    figure
    
    %plot(cam(:,3));
    hold on
    %plot(gt(:,3),'k--');
    %axis([0 inf 220 280]);
    s=VAR(sampl,inde.c(3)).^.5;
    x=MEAN(sampl,inde.c(3));
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    plot([t(1) t(end)],[1 1]*gt(3),'--k')
    set(gca,'Layer','Top')
    axis([0 inf 100 500]);
    ylabel('$c_x$ (px)')
    %     matlab2tikz('./../paper/Oleaf/fig/real_cent_x.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)
    %
    figure
    % plot(cam(:,4));
    % hold on
    % plot(gt(:,4),'k--');
    % axis([0 inf 300 350]);
    hold on
    s=VAR(sampl,inde.c(4)).^.5;
    x=MEAN(sampl,inde.c(4));
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    plot([t(1) t(end)],[1 1]*gt(4),'--k')
    set(gca,'Layer','Top')
    axis([0 inf 100 500]);
    ylabel('$c_y$ (px)')
    %     matlab2tikz('./../paper/Oleaf/fig/real_cent_y.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)
    figure
    %plot(cam(:,5));
    hold on
    s=VAR(sampl,inde.c(5)).^.5;
    x=MEAN(sampl,inde.c(5));
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    % plot([t(1) t(end)],[1 1]*x(end),'--k')
    set(gca,'Layer','Top')
    axis([0 inf -5 5]);
    ylabel('$k_1$')
    %     matlab2tikz('./../paper/Oleaf/fig/real_rad_1.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)
    figure
    
    %plot(cam(:,6));
    hold on
    s=VAR(sampl,inde.c(6)).^.5;
    x=MEAN(sampl,inde.c(6));
    fill([t(:); flipud(t(:))],[x(:)+1.96*s(:); flipud(x(:)-1.96*s(:))],1, ...
        'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(t,x,'-k')
    % plot([t(1) t(end)],[1 1]*x(end),'--k')
    set(gca,'Layer','Top')
    axis([0 inf -5 5]);
    ylabel('$k_2$')
    %     matlab2tikz('./../paper/Oleaf/fig/real_rad_2.tex', ...
    %         'height','\figureheight', ...
    %         'width','\figurewidth', ...
    %         'parseStrings', false, ...
    %         'checkForUpdates', false)


%%
if visualize
    figure
    plot3(0,0,0,'o');
    hold on
    axis equal
    for i=1:size(qte,1)
        M=quat2rmat(qte(i,:));
        vec=M*[0 0 1]';
        %plot3([0,vec(1)],[0,vec(2)],[0,vec(3)],'o-');
        plot3(vec(1),vec(2),vec(3),'kx');
        hold on
        
    end
    axis equal
    %%
    figure
    plot3(points(1,:),points(2,:),points(3,:),'xr');
    hold on
    pos=men(:,inde.p)';
    plot3(pos(1,:),pos(2,:),pos(3,:),'.k');
    
    plot3(pos(1,end),pos(2,end),pos(3,end),'oc','LineWidth',5);
    axis equal
    %%
    %[cameraParams,cameraMatrix,distCoeffs] = calibrateCameraMatlab(datapath_calib,5);
    %%
    
    cam=men(:,inde.c);
    v=ones(size(cam,1),1);
    gt=[v*cameraMatrix(1,1) v*cameraMatrix(2,2) v*cameraMatrix(1,3) v*cameraMatrix(2,3)];
    figure
    
    subplot(611)
    plot(cam(:,1));
    hold on
    plot(gt(:,1),'k--');
    axis([0 inf 400 1000]);
    title('Focal length (x) in pixels')
    subplot(612)
    plot(cam(:,2));
    hold on
    plot(gt(:,2),'k--');
    axis([0 inf 400 1000]);
    title('Focal length (y) in pixels')
    subplot(613)
    
    plot(cam(:,3));
    hold on
    plot(gt(:,3),'k--');
    axis([0 inf 100 500]);
    title('central point (x) in pixels')
    subplot(614)
    plot(cam(:,4));
    hold on
    plot(gt(:,4),'k--');
    axis([0 inf 100 500]);
    title('central point (y) in pixels')
    
    subplot(615)
    plot(cam(:,5));
    title('first radial distortion parameter')
    
    subplot(616)
    plot(cam(:,6));
    title('second radial distortion parameter')
end
if visualize
    [data,frames] = read_camera_data([folder datapath_calib],true);
end


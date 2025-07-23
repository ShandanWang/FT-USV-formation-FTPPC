function [] = AirplanePlot(varargin)

    switch(length(varargin))
        case 0
            center    = [0,0];
            direction = 0;
            axes      = gca;
            K         = 1;
        case 1
            center    = [0,0];
            direction = varargin{1};
            axes      = gca;
            K         = 1;
        case 2
            center    = varargin{1};
            direction = varargin{2};
            axes      = gca;
            K         = 1;
        case 3 
            center    = varargin{1};
            direction = varargin{2};
            axes      = gca;
            K         = varargin{3};
        case 4
            center    = varargin{2};
            direction = varargin{3};
            axes      = varargin{1};
            K         = varargin{4};
        otherwise
            error('parameters error');
    end
    if ~ishandle(axes)
        warning('Axes is not a handle; use the gca handle instead.');
        axes = gca;
    end
    if ~isnumeric(center)
        error('center is not value');
    elseif length(center)~=2
        error('center should be 1×2or2×1');
    end
    if ~isnumeric(direction)
        error('direction is not value');
    elseif length(direction)~=1
        error('direction should be value');
    end
    if ~isnumeric(K)
        error('K is not value');
    elseif length(K)~=1
        error('K should be value');
    end
%% =====================================%

    hold(axes,'on');
    if checkAB()
        filename = 'D:/2022/simulation/usv.jpg';
        [A,B] = Getlines_FromPic(filename);
        save('usv.mat','A','B')
    else
        load('C:\matlab\RFI_position\plane_data.mat','A');
        load('C:\matlab\RFI_position\plane_data.mat','B');
    end
    [A_r] = ZoomandRotate(A,direction,K);
    A_r = A_r + diag(center)*ones(size(A_r));
    multiplot(A_r,B,axes);
end

function [ flag ] = checkAB( )
    flag = 0;
    str=cd;
    files = dir(str);
    m = length(files);
    for i =1:1:m
        if strcmp(files(i).name,'usv.jpg')
            flag=1;
        end
    end
end

function [A] = ZoomandRotate(A,theta,K)
    R = [cos(theta), -sin(theta);
         sin(theta), cos(theta)];
    A=R*A*K;
end

function [A,B] = Getlines_FromPic(filename)
    frame = imread(filename);
    F_Gray=rgb2gray(frame);
    F_R=edge(F_Gray,'sobel');
    [Sx,Sy] = size(F_R);
    idx_t = int32(find(F_R));
    idx_X = idivide(idx_t,Sx)-1;
    idx_X = idx_X-mean([max(idx_X),min(idx_X)]);
    idx_Y = Sy-mod(idx_t,Sx)-1;
    idx_Y = idx_Y-mean([max(idx_Y),min(idx_Y)]);
    idx_X = double(idx_X);
    idx_Y = double(idx_Y);
    distM = max(idx_X.^2+idx_Y.^2);
    distM = sqrt(double(distM));
    idx_X = idx_X/distM/1.5;
    idx_Y = idx_Y/distM;
    [A,B] = mysort(idx_X,idx_Y,3/distM);
end

function [line_mat,Len_mat] = mysort(X,Y,th)
    Len = length(X);
    line_mat = zeros(2,Len);
    Len_mat = 1:Len;
    times = 1;
    remains_set = 2:Len;
    temp_P = [X(1),Y(1)];    
    dist_arr = zeros(1,Len); 
    th = th^2*2;
    while times<Len
        Min_dist_temp = Inf;                         
        for i = remains_set                          
            dist_temp = sum(([X(i),Y(i)]-temp_P).^2);
            if(dist_temp<Min_dist_temp)
                Min_dist_temp = dist_temp;
                Min_dist_P = i;
            end
        end
        line_mat(1,times+1)=X(Min_dist_P);
        line_mat(2,times+1)=Y(Min_dist_P);
        temp_P = [X(Min_dist_P),Y(Min_dist_P)]; 
        dist_arr(times+1) = Min_dist_temp;
        times = times+1;
        remains_set(remains_set==Min_dist_P)=[];
    end
    Len_mat(dist_arr<=th)=[];
end

function [] = multiplot(A,B,axes)
    if(nargin == 2)
        h=figure();
        axes = nexttile;
    end
    StartIdx = 2;
    for i = 1:length(B)
        temp_data_set = A(:,StartIdx:B(i)-1);
        plot(axes,temp_data_set(1,:),temp_data_set(2,:),'Color','b','LineWidth',2);
        StartIdx = B(i)+1;
        hold on;
    end
end

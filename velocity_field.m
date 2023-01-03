function out = velocity_field(imref,imdef,method,alpha,iterations,lambda,illumination,filter)
%% Velocity field calculator from Optical Flow methods
% Parameters:
%   imref   -   Reference image (2d double matrix)
%   imdef   -   Deformed image  (2d double matrix)
%   method  -   "fast" (Horn-Schunck alg.) or "slow" (Liu-Shen alg.)
%   alpha   -   a parameter that reflects the influence of the smoothness term.
%   iterations - (Horn-Shunck) number of iterations for velocity calculation
% 
% Output:
%   out     -   structure containing:
%           im1     -   the reference image
%           im2     -   the deformed image
%           u       -   u velocity field
%           v       -   v velocity field
%           U       -   sqrt(u^2+v^2) field
%           vor     -   vorticity field (not great)
%% Default input arguments
if nargin<3
    method = 'slow';
end
if nargin<4
    alpha = 20;
end
if nargin<5
    iterations = 2e2;
end
%% Parse method
switch method
    case 'fast'
        [u,v] = horn_schunck(imref,imdef);
    case 'slow'
        [u,v] = liu_shen(imref,imdef);
end
%% Conversions
% calculate the velocity magnitude
u_mag   = sqrt(u.^2+v.^2);
u_max   = max(max(u_mag));
u_mag   = u_mag/u_max;

% calculate vorticity
vor     = vorticity(u,v);
vor_max = max(max(abs(vor)));
vor     = vor/vor_max;


% calculate the 2nd invariant
Q       = invariant2_factor(u,v,1,1);

%% Output structure
out     = struct;
out.im1 = imref;
out.im2 = imdef;
out.u   = u;
out.v   = v;
out.U   = u_mag;
out.vor = vor;
%% Functions
    % Estimators
    function [u,v]  = horn_schunck(Im1,Im2)
        %% https://doi.org/10.1016/0004-3702(81)90024-2
        [u,v] = HS(Im1,Im2,alpha,iterations);
        return
        %% playground
        smoothness = [10 20 30 40];
        ui = zeros([size(Im1) numel(smoothness)]); vi = ui;
        U = @(u,v) sqrt(u.^2+v.^2);

        for i=1:numel(smoothness)
            tic
            [u,v] = HS(Im1,Im2,smoothness(i),5e2);
            t(i) = toc;
            ui(:,:,i) = u;
            vi(:,:,i) = v;
        end
        uu = U(ui,vi);
        imagesc(uu(:,:,1))
    end
    function [u,v]  = liu_shen(Im1,Im2)
        %% https://doi.org/10.1017/S0022112008003273

        %% Selete region of interest, "0" for processing the whole image, "1" for processing a selected region
        index_region=0;
        Im0 = Im1;
        Im1=double(Im1);
        Im2=double(Im2);

        if (index_region == 1)
            imagesc(uint8(Im1));
            colormap(gray);
            axis image;

            xy=ginput(2);
            x1=floor(min(xy(:,1)));
            x2=floor(max(xy(:,1)));
            y1=floor(min(xy(:,2)));
            y2=floor(max(xy(:,2)));
            I1=double(Im1(y1:y2,x1:x2));
            I2=double(Im2(y1:y2,x1:x2));
        elseif (index_region == 0)
            I1=Im1;
            I2=Im2;
        end

        I1_original=I1;
        I2_original=I2;


        % Set the Parameters for Optical Flow Computation

        % Set the lagrange multipleirs in optical computation
        lambda_1 = alpha;  % the Horn_schunck estimator for initial field
        lambda_2 = lambda; % the Liu-Shen estimator for refined estimation

        %% Number of iterations in the coarse-to-fine iterative process from
        % initial estimation, "1" means 1 iteration
        no_iteration=1; % fixed

        % Initial coarse field estimation in the coarse-to-fine iterative process,
        % scale_im is a scale factor for down-sizing of images
        scale_im=1;
        % For Image Pre-Processing

        % For local illumination intensity adjustment, To bypass it, set size_average = 0
        size_average= illumination; % in pixels

        % Gaussian filter size for removing random noise in images
        size_filter= filter; % in pixels

        % correcting the global and local intensity change in images
        [m1,n1]=size(I1);
        window_shifting=[1;n1;1;m1]; % [x1,x2,y1,y2] deines a rectangular window for global intensity correction
        [I1,I2]=correction_illumination(I1,I2,window_shifting,size_average);


        % cleaning the left and upper edges since some data near the edges are corrupted due to interperlation
        edge_width=1; % in pixels


        % pre-processing for reducing random noise,
        % and downsampling images if displacements are large
        [I1,I2] = pre_processing_a(I1,I2,scale_im,size_filter);

        I_region1=I1;
        I_region2=I2;


        % initial correlation calculation for a coarse-grained velocity field (ux0,uy0)
        % ux is the velocity (pixels/unit time) in the image x-coordinate (from the left-up corner to right)
        % uy is the velocity (pixels/unit time) in the image y-coordinate (from the left-up corner to bottom)


        % run FFT cross-correlation algorithm
        Im1=I1;
        Im2=I2;

        pivPar.iaSizeX = [64 16 8];     % size of interrogation area in X
        pivPar.iaStepX = [32 8 4];       % grid spacing of velocity vectors in X

        pivPar.ccMethod = 'fft';
        warning off
        [pivData1] = pivAnalyzeImagePair(Im1,Im2,pivPar);
        warning on

        ux0=pivData1.U;
        uy0=pivData1.V;


        % re-size the initial velocity dield (u0, v0)
        [n0,m0]=size(ux0);
        [n1,m1]=size(Im1);

        scale=round((n1*m1/(n0*m0))^0.5);

        ux0=imresize(ux0,scale);
        uy0=imresize(uy0,scale);


        % generate the shifted image from Im1 based on the initial coarse-grained velocity field (ux0, uy0),
        % and then calculate velocity difference for iterative correction


        % estimate the displacement vector and make correction in iterations

        ux=ux0;
        uy=uy0;

        k=1;
        while k<=no_iteration
            [Im1_shift,uxI,uyI]=shift_image_fun_refine_1(ux,uy,Im1,Im2);

            I1=double(Im1_shift);

            I2=double(Im2);

            % calculation of correction of the optical flow
            [dux,duy]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);

            % refined optical flow
            ux_corr=uxI+dux;
            uy_corr=uyI+duy;


            k=k+1;
        end

        % refined velocity field
        ux = ux_corr;    %%%%%
        uy = uy_corr;    %%%%%



        % clean up the edges
        ux(:,1:edge_width)=ux(:,(edge_width+1):(2*edge_width));
        uy(:,1:edge_width)=uy(:,(edge_width+1):(2*edge_width));

        ux(1:edge_width,:)=ux((edge_width+1):(2*edge_width),:);
        uy(1:edge_width,:)=uy((edge_width+1):(2*edge_width),:);

        u = ux;
        v = uy;
    end
    function [u, v] = HS(im1, im2, alpha, ite, uInitial, vInitial, displayFlow, displayImg)
        % Horn-Schunck optical flow method
        % Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
        % 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337
        %
        % Usage:
        % [u, v] = HS(im1, im2, alpha, ite, uInitial, vInitial, displayFlow)
        % For an example, run this file from the menu Debug->Run or press (F5)
        %
        % -im1,im2 : two subsequent frames or images.
        % -alpha : a parameter that reflects the influence of the smoothness term.
        % -ite : number of iterations.
        % -uInitial, vInitial : initial values for the flow. If available, the
        % flow would converge faster and hence would need less iterations ; default is zero.
        % -displayFlow : 1 for display, 0 for no display ; default is 1.
        % -displayImg : specify the image on which the flow would appear ( use an
        % empty matrix "[]" for no image. )
        %
        % Author: Mohd Kharbat at Cranfield Defence and Security
        % mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
        % Published under a Creative Commons Attribution-Non-Commercial-Share Alike
        % 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
        %
        % October 2008
        % Rev: Jan 2009

        %% Default parameters
        if nargin<1 || nargin<2
            im1=imread('yos9.tif');
            im2=imread('yos10.tif');
        end
        if nargin<3
            alpha=1;
        end
        if nargin<4
            ite=100;
        end
        if nargin<5 || nargin<6
            uInitial = zeros(size(im1(:,:,1)));
            vInitial = zeros(size(im2(:,:,1)));
        elseif size(uInitial,1) ==0 || size(vInitial,1)==0
            uInitial = zeros(size(im1(:,:,1)));
            vInitial = zeros(size(im2(:,:,1)));
        end
        if nargin<7
            displayFlow=0;
        end
        if nargin<8
            displayImg=im1;
        end

        %% Convert images to grayscale
        if size(size(im1),2)==3
            im1=rgb2gray(im1);
        end
        if size(size(im2),2)==3
            im2=rgb2gray(im2);
        end
        im1=double(im1);
        im2=double(im2);

        im1=smoothImg(im1,1);
        im2=smoothImg(im2,1);

        tic;

        %%
        % Set initial value for the flow vectors
        u = uInitial;
        v = vInitial;

        % Estimate spatiotemporal derivatives
        [fx, fy, ft] = computeDerivatives(im1, im2);
        % [fxb, fyb, ftb] = computeDerivatives(im1, im2,'barron');
        % [fxd, fyd, ftd] = computeDerivatives(im1, im2,'diff');

        % Averaging kernel
        kernel_1=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];

        % Iterations
        % diffs = zeros(1,ite);
        for i=1:ite
            % Compute local averages of the flow vectors
            uAvg=conv2(u,kernel_1,'same');
            vAvg=conv2(v,kernel_1,'same');
            %     u0 = u;
            %     v0 = v;
            % Compute flow vectors constrained by its local average and the optical flow constraints
            u= uAvg - ( fx .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
            v= vAvg - ( fy .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
            %     diffs(i) = mean2(sqrt(u.^2+v.^2)-sqrt(u0.^2+v0.^2));
            %     title(sprintf('iterations = %i',i))
            %     pause(.01)
        end

        % plot(1:ite,diffs(1)-diffs)
        u(isnan(u))=0;
        v(isnan(v))=0;

        %% Plotting
        if displayFlow==1
            plotFlow(u, v, displayImg, 5, 5);
        end
    end
    % Calculators
    function [fx, fy, ft]   = computeDerivatives(im1, im2,method)
        if nargin<3
            method = 'hs';
        end

        if size(im2,1)==0
            im2=zeros(size(im1));
        end

        switch method
            case 'hs'
                % Horn-Schunck original method
                fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
                fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
                ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
            case 'barron'
                % derivatives as in Barron
                fx= conv2(im1,(1/12)*[-1 8 0 -8 1],'same');
                fy= conv2(im1,(1/12)*[-1 8 0 -8 1]','same');
                ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
                fx=-fx;fy=-fy;
            case 'diff'
                % An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
                fx = conv2(im1,[1 -1],'same');
                fy = conv2(im1,[1; -1],'same');
                ft= im2-im1;
        end
    end
    
    function [QQ]           = invariant2_factor(Vx, Vy, factor_x, factor_y)
        % factor_x: converting factor from pixel to m (m/pixel) in x
        % factor_y: converting factor from pixel to m (m/pixel) in y

        % Vx = imfilter(Vx, [1 1 1 1 1]'*[1 1 1 1 1]/25,'symmetric');
        % Vy = imfilter(Vy, [1 1 1 1 1]'*[1,1 1 1,1]/25,'symmetric');

        dx=1;
        D = [0, -1, 0; 0,0,0; 0,1,0]/2; %%% partial derivative
        Vx_x = imfilter(Vx, D'/dx, 'symmetric',  'same')/factor_x;
        Vx_y = imfilter(Vx, D/dx, 'symmetric',  'same')/factor_y;

        Vy_x = imfilter(Vy, D'/dx, 'symmetric',  'same')/factor_x;
        Vy_y = imfilter(Vy, D/dx, 'symmetric',  'same')/factor_y;


        [M,N]=size(Vx);

        for m=1:M
            for n=1:N
                uu(1,1)=Vx_x(m,n);
                uu(1,2)=Vx_y(m,n);
                uu(2,1)=Vy_x(m,n);
                uu(2,2)=Vy_y(m,n);

                for i=1:2
                    for j=1:2
                        S(i,j)=(uu(i,j)+uu(j,i))/2;
                        Qq(i,j)=(uu(i,j)-uu(j,i))/2;
                    end
                end

                QQ(m,n)=(trace(Qq*Qq')-trace(S*S'))/2;


                %         lambda1(m,n)=max(d);
                %         lambda2(m,n)=min(d);
            end
        end

    end
    % Plotting
    function plotFlow(u, v, imgOriginal, rSize, scale)
            % Creates a quiver plot that displays the optical flow vectors on the
            % original first frame (if provided). See the MATLAB Function Reference for
            % "quiver" for more info.
            %
            % Usage:
            % plotFlow(u, v, imgOriginal, rSize, scale)
            %
            % u and v are the horizontal and vertical optical flow vectors,
            % respectively. imgOriginal, if supplied, is the first frame on which the
            % flow vectors would be plotted. use an empty matrix '[]' for no image.
            % rSize is the size of the region in which one vector is visible. scale
            % over-rules the auto scaling.
            %
            % Author: Mohd Kharbat at Cranfield Defence and Security
            % mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
            % Published under a Creative Commons Attribution-Non-Commercial-Share Alike
            % 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
            %
            % October 2008
            % Rev: Jan 2009

            figure();

            if nargin>2
                if sum(sum(imgOriginal))~=0
                    imshow(imgOriginal,[0 255]);
                    hold on;
                end
            end
            if nargin<4
                rSize=5;
            end
            if nargin<5
                scale=3;
            end

            % Enhance the quiver plot visually by showing one vector per region
            for i=1:size(u,1)
                for j=1:size(u,2)
                    if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
                        u(i,j)=0;
                        v(i,j)=0;
                    end
                end
            end
            quiver(u, v, scale, 'color', 'b', 'linewidth', 2);
            set(gca,'YDir','reverse');
    end
end
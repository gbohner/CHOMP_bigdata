function update_visualize( y,H, W, opt,cur_pc,varargin)
%UPDATE_VISUALIZE Summary of this function goes here
%   Detailed explanation goes here



if nargin >5
  show_numbers = varargin{1};
else
  show_numbers = 0;
end

if nargin >6
  show_specific = varargin{2};
else
  show_specific = [];
end

NSS = opt.NSS;
KS = opt.KS;
Nmaps = size(W,3);
isfirst = zeros(1,Nmaps);
m = size(W,1);
d = (m-1)/2;
xs  = repmat(-d:d, m, 1);
ys  = xs';
rs2 = (xs.^2+ys.^2);

% Sets the opacity of individual red dots according to their position in the order
alphaspace = logspace(1,0,numel(H))/10; %TODO make that it is according to the score instead of position in order (needs more input to the func tho)

    sign_center = -squeeze(sign(W(d,d,:)));
    sign_center(:) = 1;
    Wi = reshape(W, m^2, Nmaps);
    nW = max(abs(Wi), [], 1);
    %             nW = sum(Wi.^2, 1).^.5;

    Wi = Wi./repmat(sign_center' .* nW, m*m,1);

	
    figure(1); 
			if cur_pc == 0, set(gcf, 'Visible', 'off'), end
			visualSS(Wi, 4, KS, [-1 1]); colormap('jet')
			%if cur_pc == 0, print([opt.output_folder filesep opt.output_file_prefix '_fig1.eps'],'-depsc2'), end

    figure(3);
						if cur_pc == 0, set(gcf, 'Visible', 'off'), end
	colormap('jet')
%         Im = y(:,:,ex);
    Im = y;
    sig = nanstd(Im(:)); mu = nanmean(Im(:)); M1= mu - 4*sig; M2= mu + 12*sig;
    imagesc(Im, [M1 M2]);
    colormap gray
    mycolor = 'rymg';



    if isempty(show_specific)
      show_specific = 1:size(H,1);
    end
    hold on
    axis image
    for i12 = show_specific%[1,5,10,20,30, 50] %1:length(col)
      [row,col,type] = ind2sub(size(y),H(i12));
      if show_numbers
        text(col, row, num2str(i12), 'Color',mycolor(mod(type-1,length(mycolor))+1),'FontSize',20,'FontWeight','bold');
      else
        %tmp1 = plot(col, row, 'or' , 'Linewidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', mycolor(mod(type-1,length(mycolor))+1), 'MarkerEdgeColor', mycolor(mod(type-1,length(mycolor))+1));
        tmp1 = scatter(col, row, 30, 'r', 'filled'); %TODO show both multiple colors and alphas
        alpha(tmp1, alphaspace(i12))
      end
%             text(col(i12), row(i12), num2str(i12), 'Color',mycolor(mod(map(i12)-1,length(mycolor))+1));
%             title(i12);
%             waitforbuttonpress; 
    end

    hold off

    drawnow
	%		if cur_pc == 0, print([opt.output_folder filesep opt.output_file_prefix '_fig3.eps'],'-depsc2'), end

      
  function visualSS( W, mag, rows,clims, varargin)
    % visual - display a basis for image patches
    %
    % W        the basis, with patches as column vectors
    % mag      magnification factor
    % cols     number of columns (x-dimension of map)
    % ysize    [optional] height of each subimage
    %

    mini=min(W(:));

    % This is the side of the window

    if ~isempty(varargin)
        A = varargin{1};    
        xsize = A(1);
        ysize = A(2);
    else
        ysize = sqrt(size(W,1));
        xsize = size(W,1)/ysize;
    end

    % Helpful quantities
    xsizem = xsize-1;
    xsizep = xsize+1;
    ysizem = ysize-1;
    ysizep = ysize+1;
    cols = ceil(size(W,2)/rows);

    % Initialization of the image
    I = mini*ones(2+ysize*rows+rows-1,2+xsize*cols+cols-1);

    for j=0:cols-1
        for i=0:rows-1        
            if j*rows+i+1>size(W,2)
                1;
                % This leaves it at background color            
            else
                % This sets the patch
                I(i*xsizep+2:i*xsizep+xsize+1, ...
                    j*ysizep+2:j*ysizep+ysize+1) = ...
                    reshape(W(:,j*rows+i+1),[xsize ysize]);
            end

        end
    end

    % Make a black border
    I(1,:) = 0;
    I(:,1) = 0;
    I(end,:) = 0;
    I(:,end) = 0;

    I = imresize(I,mag);

    imagesc(I, clims);

    axis off
    colormap('gray')
    % colorbar
    % iptsetpref('ImshowBorder','tight'); 
    truesize;  
    drawnow
  end
      
end


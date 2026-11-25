function [im, s] = readIm(folder, file, widMax) % read image from image file

	%	Read file and skip ineligible images
	file = folder + filesep + file; % file name with full path
	try % some files are unreadable
		im = imread(file); % image: ys x xs x cs, uint8
	catch % error when reading
		im = []; s = []; % return empty image and size
		return
	end
	s = size(im); % original image size
	if numel(s) < 3 % greyscale file
		im = []; % return empty image
		s = [s, 1]; % return size
		return 
	end

	%	Crop images that are too big and crop to square
	if nargin > 2 % lim exists
		l = widMax; % upper limit on size (pix): 1 x 1
		l = min(s, l); % limits consistent with current size: 1 x 3
		l = min(l(1: 2)); % smaller of width and height: 1 x 1
		l = [l, l]; % make it square: 1 x 2
		w = centerCropWindow2d(s, l); % cropping window
		im = imcrop(im, w); % crop image
	end

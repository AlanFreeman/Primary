function [name, val, vals] = getIm(m, task) % get images to be used as stimuli

	% Obtain images
	load(m.(task).file, 'd'); % table containing image, one row per image
	d = removevars(d, 'wave'); % d.wave is large and not required
	i = m.(task).index; % indices of images to retain
	if ~ isempty(i) % use all images if i is empty
		d = d(i, :); % select required images
	end
	d = d(d.valid == true, :); % keep only valid files
	is = size(d, 1); % number of images

	switch m.p.stimT % drifting or pulsed image?
		case 'drift' % drifting image
			%{
			name = {'file', 'dir', 'win'}; % names of stimulus parameters
			win = m.(task).win;
			val = {d.image, dir, win}; % values of all stimulus parameters
			vals = [is, length(dir), length(win)]; % number of stimulus parameters
			%}
			name = {'index', 'dir'}; % names of stimulus parameters
			dir = m.(task).dir; % direction (deg): ds x 1
			val = {d.image, m.(task).dir}; % values of all stimulus parameters
			vals = [is, length(dir)]; % number of stimulus parameters
		case 'pulse' % pulsed image

			%	Crop and resize each image
			widS = 101; % final size of image, width and heights (pix)
			im = zeros(widS, widS, 3, is, 'uint8'); % storage: ys x xs x cs x is
			for i = 1: is % loop over images

				%	Read table
				dC = d(i, :); % current row of table
				imC = dC.image.im; % current image (sRGB-units): ys x xs x cs, uint8
				s = size(imC); wid = s(1); % width and height (pix)
				waveDom = dC.waveDom; % dominate wavelength in image (pix/cycle)

				%	Crop the image to match receptive field size, and to standardise
				cycIm = wid / waveDom; % number of cycles in image
				cycField = m.p.freqSPref * m.p.wid; % number of cycles in vis. field
				r = cycField / cycIm; % cycle ratio
				if r > 1 % can't crop because image is too small
					continue
				end 
				wid = floor(r * wid); % reduced width to match receptive field
				win = centerCropWindow2d(s, [wid, wid]); % cropping window
				imC = imcrop(imC, win); % crop image to match receptive field sizes
				imC = imresize(imC, [widS, widS]); % resize to standard
				im(:, :, :, i) = uint8(imC); % store: ys x xs x cs x is, uint8

			end

			% Reformat image to standard array, remove small images, store
			im = permute(im, [2, 1, 3, 4]); % map: xs x ys x cs x is, uint8
			im = flip(im, 2); % map: xs x ys x cs x is, uint8
			name = {'image'}; % place holder
			val = im; % images: xs x ys x cs x is, uint8
			vals = size(im, 4); % number of images

	end
	
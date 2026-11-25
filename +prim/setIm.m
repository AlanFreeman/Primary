function m = setIm(m, val, sub) % set up the image stimulus

%	Inputs:
%		m = metadata
%		name = names of specified stimulus parameters: cell
%		val = values for the specified stimulus parameters:
%			cell array with one member for each stimulus parameter
%		sub = subscripts of current set of values: cell

	import prim.see % find function
	switch m.p.stimT % drifting or pulse image?
		case 'drift' % construct a movie by panning across the current image

			%	Obtain image and convert it to model format
			im = val{1}; subC = sub{1}; % values and subscript
			im = im(subC).im; % current image (uint8): ys x xs x cs
			%	see(m, image = im, format = 'imshow') % debug
			im = permute(im, [2, 1, 3]); % swap x and y: xs x ys x cs, uint8
			im = flip(im, 2); % invert y: xs x ys x cs, uint8
			im = rgb2lin(im); % linearise gamma-corrected values: xs x ys x cs, uint8
			%	see(m, image = im) % debug

			%	Obtain movie parameters
			valC = val{2}; subC = sub{2}; % values and subscript
			dir = pi * valC(subC) / 180; % direction in which to pan (rad)
			%{
			valC = val{3}; subC = sub{3}; % values and subcript
			win = valC(subC); % fraction of image to be shown in each movie frame
			%}
			win = .75; % fraction of image to be shown in each movie frame

			% Spatially attenuate the image with the optical point spread function
			sig = (1 / sqrt(2)) * m.p.radOpt; % convert PSF radius to st. dev. (deg)
			wI = size(im, 1); % image width (px)
			wF = round(win * wI); % frame width (px)
			sig = (sig / m.p.wid) * wF; % PSF standard deviation (pixels)
			im = imgaussfilt(im, sig); % cross-corr. image, PSF: xs x ys x cs, uint8
			%	see(m, image = im) % debug

			%	Find frame centrepoint at end of travel
			t = tan(dir); % tan of direction
			if abs(t) <= 1 % endpoint is on vertical boundary
				xE = [1, t]; % endpoint (travel fraction): 1 x 2
			else % endpoint is on horizontal boundary
				xE = [1 / t, 1]; % endpoint (travel fraction): 1 x 2
			end
			if sin(pi / 4 + dir) < 0, xE = - xE; end % endpoint is below neg. diagonal

			%	Calculate frame centre during movie
			ts = m.p.ts; % number of frames
			x = linspace(-1, 1, ts)'; % 1D frame centre (travel fraction): ts x 1
			x = xE .* x; % 2D frame centre (travel fraction): ts x 2

			%	Create movie
			wT = wI - wF; % travel distance across image (px)
			x = .5 * wT * x; % frame centre relative to image centre (px): ts x 2
			x = .5 * wT + x; % frame origin relative to image origin (left bottom)
			j = 1: wF; % indices of pixels in frame: 1 x ps
			imF = zeros(wF, wF, 3, ts, 'uint8'); % allocate storage: xs x ys x cs x ts
			for i = 1: ts % loop over frames
				xC = round(x(i, :)); % current frame origin (px): 1 x 2
				imC = im(xC(1) + j, xC(2) + j, :); % current frame (uint8): xs x ys x cs
				imF(:, :, :, i) = imC; % crop image to frame				
			end
			%	see(m, movie = imF) % debug
			m.p.image = imF; % store image: xs x ys x cs x ts, uint8

		case 'pulse' % pulsed image *** convert iVal to sub ***
			m.p.image = val(:, :, :, iVal); % current image (uint8): xs x ys x cs
	end

function [mag, t] = filtIm(im) % create Gabor filter bank and filter image

	% Define wavelengths and directions
	[ys, xs, ~] = size(im); % image size: ys x xs x cs
	wMax = .5 * hypot(ys, xs); % maximum wavelength (pix)
	switch 'lin' % wavelength spacing
		case 'lin' % linear
			ws = 10; % number of wavelengths
			w = linspace(0, wMax, ws + 1); % wavelength (pix)
			w(1) = []; % exclude zero
		case 'log' % logarithmic
			wMin = 2; % minimum wavelength (pix)
			ws = ceil(log2(wMax / wMin)); % number of wavelengths
			w = 2 .^ (- ws + 1: 0) * wMax; % wavelengths (pix): 1 x ws
	end
	d = linspace(-90, 90, 10 + 1); % stim. d'n, closed interval (deg)
	d = d(1: end - 1); % stimulus direction, open interval (deg)

	%	Create table of wavelengths and directions, and Gabor filter bank
	[wT, dT] = ndgrid(w, d); % grid of w, d: ws x ds
	t = [wT(:), dT(:)]; % table of w, d: gs x 2, gs = ws * ds
	b = gabor(w, d); % Gabor filter bank: 1 x gs, gabor object
	
	%	Calculate factor to normalise kernels
	n = arrayfun(@(x) abs(x.SpatialKernel), b, 'uniformOutput', false);
		% absolute value of kernels: 1 x gs, cell
	n = cellfun(@(x) sum(x(:)), n); % sum of absolute values: 1 x gs
	n = shiftdim(n, -1); % prepare for division: 1 x 1 x gs

	% Apply Gabor filters and find indices of maximum magnitude
	im = im(:, :, 1); % Red component only: most salient component in macaques
	mag = imgaborfilt(im, b); % Gabor magnitude: ys x xs x gs
	mag = mag ./ n; % normalise magnitudes: ys x xs x gs

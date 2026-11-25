function [mag, t] = filtImGpu(im) % filter image using GPU

%	Takes longer than version without GPU: not sure why

	% Define wavelengths and directions
	im = im(:, :, 1); % Red component only: most salient component in macaques
	im = im2double(im); % convert to double for imfilter: ys x xs
	[ys, xs, ~] = size(im); % image size: ys x xs x cs
	im = gpuArray(im); % shift to GPU
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
	ds = 10; % number of directions
	d = linspace(-90, 90, ds + 1); % direction, closed interval (deg)
	d = d(1: end - 1); % direction, open interval (deg)

	%	Create table of wavelengths and directions, and Gabor filter bank
	[wT, dT] = ndgrid(w, d); % grid of w, d: ws x ds
	t = [wT(:), dT(:)]; % table of w, d: gs x 2, gs = ws * ds
	b = gabor(w, d); % Gabor filter bank: 1 x gs, gabor object
	%	b = reshape(b, [ws, ds]); % Gabor filter bank: ws x ds
	gs = length(b); % number of Gabors

	%	Loop over Gabors
	mag = zeros(ys, xs, gs); % allocate storage
	for i = 1: gs % loop over Gabors
		%	w = t(i, 1); d = t(i, 2); % current wavelength and direction
		g = b(i); % current Gabor object: 1 x 1, struct
		k = g.SpatialKernel; % Gabor: complex double
		imReal = imfilter(im, real(k), 'replicate'); % filt. image, real: ys x xs
		imImag = imfilter(im, imag(k), 'replicate'); % filt. image, imag.: ys x xs
		magC = abs(complex(imReal, imImag)); % magnitude of filtered image: ys x xs
		n = abs(k); % normalising factor: ys x xs
		n = sum(n(:)); % sum of absolute values: 1 x 1
		mag(:, :, i) = magC ./ n; % normalise magnitudes: ys x xs x gs
	end
	mag = gather(mag); % shift to the CPU

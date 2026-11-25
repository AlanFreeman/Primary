function d = transIm(m) % transduce image: interpolate, convert to contrast

	%	Initialise
	import prim.see prim.rgb2cont % find functions
	im = m.p.image; % image: xs x ys x cs x ts, uint8
	%	see(m, image = im) % debug
	[xs, ~, ~, ts] = size(im); % image dimensions
	loc = m.p.cone.loc; % cone location (deg): ls x 2
	ls = size(loc, 1); % number of cones	
	type = m.p.cone.type; % cone type: 1, 2, 3 for L, M, S

	%	Transduce
	switch m.p.stimT % drifting or pulsed image?
		case 'drift' % drifting image

			%	Convert the image to cone contrast
			im = rgb2cont(im, m, 'linear');
				% convert RGB to cone contrast: xs x ys x cs x ts, double

			%	Interpolate at cone locations
			x = .5 * m.p.wid * linspace(-1, 1, xs); % x locations (deg): 1 x ps
			locX = loc(:, 1); locY = loc(:, 2); % cone location (deg): 1 x ls
			im = interpn(x, x, im, locX, locY); % interpolate: ls x 1 x cs x ts
		
			%	Select cone-specific stimulus
			d = zeros(ls, ts); % allocate storage: ls x ts
			for i = 1: ls % loop over cones
				j = type(i); % cone type
				d(i, :) = im(i, :, j, :); % select appropriate cone contrast: ls x ts
			end
		
			%	Transform the drive signal to the temporal frequency domain
			w = hann(ts, 'periodic'); % Hann window
			d = w' .* d; % make the drive periodic
			d = fft(d, [], 2); % Fourier transform (contrast-units): ls x ts

		case 'pulse' % pulsed image *** fix ***

			%	Convert the image to cone contrast
			im = rgb2cont(im, m); % convert sRGB to cone cont.: xs x ys x cs, double
			%	image(.5 * (1 + flipud(permute(im, [2, 1, 3]))))
		
			% Spatially attenuate the image with the optical point spread function
			sig = (1 / sqrt(2)) * m.p.radOpt; % convert PSF radius to st. dev. (deg)
			sig = (sig / m.p.wid) * xs; % point spread function radius (pixels)
			im = imgaussfilt(im, sig); % cross-corr. image, PSF: xs x ys x cs
		
			%	Interpolate at cone locations
			x = .5 * m.p.wid * linspace(-1, 1, xs); % x locations (deg): 1 x xs
			locX = loc(:, 1); locY = loc(:, 2); % cone location (deg): 1 x ls
			v = string(version('-release')); % Matlab version number
			if v < "2021a" % old versions can't process cs > 1
				d = zeros(length(locX), 1, 3); % drive: ls x 1 x cs
				for i = 1: 3 % loop over colours
					d(:, 1, i) = interpn(x, x, im(:, :, i), locX, locY);
						% interpolate: ls x 1 x cs
				end
			else
				d = interpn(x, x, im, locX, locY); % interpolate: ls x 1 x cs
			end

			%	Add temporal component and select appropriate contrast
			dC = d .* pulse(m); % multiply by pulse transform (mV): ls x fs x cs
			d = zeros(ls, fs); % allocate storage: ls x fs
			for i = 1: ls % loop over cones
				j = type(i); % cone type
				d(i, :) = dC(i, :, j); % select appropriate cone contrast: ls x fs
			end

	end
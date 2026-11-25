function im = rgb2cont(im, m, linearity) % convert sRGB to cone contrast

	% Convert RGB to XYZ
	import prim.see % find function
	switch string(linearity)
		case "linear" % assume that gamma correction has already been removed
			im = rgb2xyz(im, colorSpace = "linear-rgb"); % convert RGB to XYZ: double
		case "nonlinear" % remove gamma correction
			im = rgb2xyz(im); % convert sRGB to XYZ: double
	end
	%	see(m, image = im) % debug

	%	Convert XYZ to LMS
	[xs, ys, cs, ts] = size(im); % size of image
	im = permute(im, [3, 1, 2, 4]); % prepare for mult.: cs x xs x ys x ts
	im = reshape(im, cs, []); % prepare for mult.: cs x ps, ps = xs * ys * ts
	im = m.p.xyz2lms * im; % convert XYZ to LMS: (cs x cs) * (cs x ps) = cs x ps
	im = reshape(im, [cs, xs, ys, ts]); % cs x xs x ys x ts
	im = permute(im, [2, 3, 1, 4]); % xs x ys x cs x ts
	%	see(m, image = im) % debug
	%	see(m, movie = im) % debug
	
	%	Convert LMS to cone contrast
	sig = (1 / sqrt(2)) * m.p.radFix ; % convert radius to standard dev. (deg)
	sig = (xs / m.p.wid) * sig; % background radius (pixels)
	b = imgaussfilt(im, sig); % local background: xs x ys x cs x ts
	im = (im - b) ./ b; % contrast: xs x ys x cs x ts, double
	%	see(m, image = im, format = "imageCont") % debug
	%	see(m, movie = im) % debug

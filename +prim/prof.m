function [d, m] = prof(d, m)

%	Create a one-dimensional profile from a receptive field

	%	Label the data file rows
	cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1]; % contrast
	stim = ["M", "equi", "L", "ach"]; % stimulus types
	for i = 1: size(d, 1) % loop over data file rows
		j = ismember(cont, d.cont(i, :), "rows"); % which stimulus type?
		d.Properties.RowNames(i) = stim(j); % assign type
	end
	d = d(["M", "L"], :); % select L and M rows: rs x vs
	d = sortrows(d, "RowNames"); % sort table by row name: L then M

	%	Read and reshape response
	r = d.resp; % response (Hz): rs x 1 x ps x cs
	rL = r(1, 1, :, 1); rM = r(2, 1, :, 2); % select corresponding colours
	r = [rL; rM]; % recombine: rs x 1 x ps
	[~, ~, ps] = size(r); % number of pixels
	xs = sqrt(ps); % number of x values in receptive field
	r = reshape(r, [2, xs, xs]); % reshape to form image: rs x xs x ys

	% Rotate the receptive field so that the preferred direction is parallel to x
	if isfield(m.prof, "angle") % set preferred direction
		a = m.prof.angle; % user-specified angle (deg)
	else
		a = 0; % default angle (deg)
	end
	if a % rotation is required
		r = shiftdim(r, 1); % reshape for rotation: xs x ys x rs
		%	prim.see(m, image = r(:, :, 1), format = "imageCont") % debug
		r = imrotate(r, - a, "nearest", "crop");
			% rotate anticlockwise by required angle
		%	prim.see(m, image = r(:, :, 1), format = "imageCont") % debug
		r = shiftdim(r, 2); % return to original shape: rs x xs x ys
	end

	%	Sum across y values to obtain profile, and normalise
	r = sum(r, 3); % sum across y: rs x xs
	n = max(abs(r), [], "all"); % largest absolute value
	r = r / n; % normalise

	%	Store
	x = .5 * m.p.wid * linspace(-1, 1, xs); % x values (deg)
	d.x = repmat(x, [2, 1]); % store x values: rs x xs
	d.resp = r; % store response: rs x xs
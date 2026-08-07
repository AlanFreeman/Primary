function [dens, x] = densFunc(m) % plot density map of L and M cones

%	Method:
%		The map is designed for comparison with cortical receptive fields.
%		Cone locations are therefore cross-correlated with a Gaussian representing
%		the convergence function between stimulus and cortex.

	% Calculate map locations
	if isfield(m.dens, 'xs') % number of x locations in map
		xs = m.dens.xs; % specified by user
	else
		xs = 101; % default
	end
	ls = xs ^ 2; % number of locations in map

	%	Create a list of map locations
	x = .5 * m.p.wid * linspace(-1, 1, xs); % x values (deg): 1 x xs
	[xG, yG] = ndgrid(x); % map locations (deg): xs x ys
	locMap = [xG(:), yG(:)]; % list of map locations (deg): ls x 2

	%	Insert cones
	loc = m.p.cone.loc; % cone locations (deg): cs x 2
	t = m.p.cone.type; % cone type: L, M, S
	dens = zeros(ls, 2); % map in list form: columns 1, 2 for L-, M-cone, resp. 
	for i = 1: 2 % L then M
		locC = loc(t == i, :); % cone locations (deg): cs x 2
		j = knnsearch(locMap, locC); % index of cone locations in map
		dens(j, i) = 1; % insert cones: ls x 2
	end
	dens = reshape(dens, [xs, xs, 2]); % reshape to map: xs x ys x 2

	%	Cross-correlate map with Gaussian convergence function
	sig = m.p.radBetaCen / sqrt(2); % convert radius to standard deviation (deg)
	sig = (sig / m.p.wid) * xs; % convergence function radius (pixels)
	dens = imgaussfilt(dens, sig); % density (cells/v.f.): xs x ys x 2

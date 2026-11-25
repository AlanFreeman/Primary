function [d, m] = crop(d, m)

% Crop locations to smaller field
%
%	Input: table can have multiple rows

	% Prepare response
	r = d.(m.z); % response: ps x ls x qs; ps and qs can be multi-dimensional
	loc = shiftdim(d.loc(1, :, :), 1); % location (deg): ls x 2
	if isfield(m.crop, "radius") % radius is specified
		rad = m.crop.radius; % radius (deg)
	else % default radius
		rad = m.p.radSur; % surround radius (deg)
	end
	rad = .5 * m.p.wid - rad; % half-width of cropped field (deg)
	if isfield(m.crop, "dim") % set location dimension
		dim = m.crop.dim; % user specified
	else % default
		dim = d.Properties.CustomProperties.RespDim; % response dimensions
		[~, dim] = ismember("loc", dim); % index of location in response
	end

	%	Crop and store
	[r, loc] = doCrop(r, loc, rad, dim);
	d.(m.z) = r; % store response: ps x ls x qs, ps and qs can be multi-dim.
	loc = shiftdim(loc, -1); % location (deg): 1 x ls x 2
	d.loc = repmat(loc, [size(d, 1), 1, 1]); % location (deg): rs x ls x 2

function [r, loc] = doCrop(r, loc, rad, dim) % crop border locations

%	Input:
%		r = response: multi-dimensional
%		loc = response locations: ls x 2
%		rad = half-width of square border in which to retain responses
%		dim = location dimension in r
%	Output:
%		r = cropped response
%		loc = retained locations

	s = size(r); % size of r
	es = prod(s(1: dim - 1)); % number of elements before dimension dim
	r = reshape(r, es, s(dim), []); % prepare for cropping
	i = abs(loc(:, 1)) <= rad & abs(loc(:, 2)) <= rad; % indices of loc. to keep
	r = r(:, i, :); % crop
	s(dim) = sum(i); % update size
	r = reshape(r, s); % restore original shape
	loc = loc(i, :); % update locations
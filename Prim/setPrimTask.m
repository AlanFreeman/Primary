function m = setPrimTask % define function handles for setPrim tasks

	h = localfunctions; % handles of functions in this file
	hs = length(h); % number of handles
	for i = 1: hs % loop over handles
		hC = h{i}; % current handle
		task = func2str(hC); % name of task
		m.(task).fun = hC; % store handle in metadata
	end

function [d, m] = countCen(d, m) % cones per ganglion cell centre

	% Obtain radius data
	[dC, m] = rad(d, m); % calculate radii
	i = dC.source == 'cenConv'; % indices for centre radii
	dC = dC(i, :); % select centre data: es x 1
	r = dC.radius; % centre radius (deg): es x 1
	es = size(dC, 1); % number of eccentricities

	% Obtain density data from input table
	ecc = d.eccDeg; % eccentricity (deg): eDs x 1, eDs = number of densities
	dens = d.densDeg; % cone density (cones/deg^2): eDs x 1
	dens = interp1(ecc, dens, dC.ecc); % interp. density data (deg^-2): es x 1
	c = pi * r .^ 2 .* dens; % area * density = cone count: es x 1
	d = repmat(d(1, :), [es, 1]); % change number of rows to es
	d.eccDeg = dC.ecc; d.count = c; % store
	d.Properties.VariableDescriptions{'count'} = 'Cones per ganglion cell centre';

function [d, m] = countCone(d, m) % cone count vs eccenticity: Packer (89)

	i = d.quad == 'temporal'; d = d(i, :); % keep only temporal quadrant
	e = d.eccDeg; % eccentricity (deg): es x 1
	dens = d.densDeg; % density (deg^-2): es x 1
	dn = 2 * pi * e .* dens; % derivative of cell count (deg^-1): es x 1
	n = cumtrapz(e, dn); % integrate to find cell count: es x 1
	d.count = n; % store cell count
	d.Properties.VariableDescriptions{'count'} = 'Cone count';

function [d, m] = countGang(d, m) % g.c. count vs fun. ecc.: Wässle (89)

	% Initialise
	i = d.type == 'func'; % functional eccentricity
	d = d(~ i, :); % keep only anatomical eccentricity
	d = sortrows(d, 'eccDeg'); % sort on eccentricity
	eA = d.eccDeg; % anatomical eccentricity (deg): es x 1
	if isfield(m.p, 'densGangCoef') % regression coefficients are available
		densA = exp(polyval(m.p.densGangCoef, log10(eA))); % dens. (deg^-2): es x 1
	else
		densA = d.densDeg; % density (deg^-2): es x 1
	end
	
	% Integrate density to find cell count
	switch 'cont' % method
		case 'cont' % continuous
			dn = 2 * pi * eA .* densA; % derivative of cell count (deg^-1): es x 1
			n = cumtrapz(eA, dn); % integrate to find cell count: es x 1
		case 'disc' % discrete
			dE = diff(eA); % eccentricity increments (deg); (es - 1) x 1
			dE = [dE(1); dE]; % replicate first point (deg): es x 1
			a = pi * (eA - .5 * dE) .^ 2; % circle areas (deg^2): es x 1
			a = [a; pi * (eA(end) + .5 * dE(end)) ^ 2]; % add final pt.: (es + 1) x 1
			dA = diff(a); % annulus area (deg^2): es x 1
			n = dA .* densA; % annulus count: es x 1
			n = cumsum(n); % cumulative count: es x 1
	end

	% Convert anatomical to functional eccentricity, and store
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	e = interp1([0; dE.ecc], [0; dE.eccFun], eA);
		% functional eccentricity (deg), extended to origin: es x 1
	d.ecc = e; % store functional eccentricity
	d.count = n; % store ganglion cell count
	i = d.eccDeg == 0 | d.eccDeg >= dE.ecc(1); % rows in conversion range
	d = d(i, :); % keep them

	%	Describe
	d = renamevars(d, 'eccDeg', 'eccAnat'); % rename for clarity
	d.Properties.VariableDescriptions{'ecc'} = 'Functional eccentricity (deg)';
	d.Properties.VariableDescriptions{'eccAnat'} = ...
		'Anatomical eccentricity (deg)';
	d.Properties.VariableDescriptions{'count'} = 'Ganglion cell count';

function [d, m] = cover(d, m) % ganglion cell coverage: multiple authors

% Input: ganglion cell density versus functional eccentricity

	dens = d.densDeg; % density (cells/deg^2)
	switch 'cen' % source of radius
		case 'cen' % ganglion cell centre radius
			dC = rad(d, m); % calculate centre radius
			i = dC.source == 'cenConv'; dC = dC(i, :); % select centre radius
			r = dC.radius; % centre radius (deg)
			dens = interp1(d.ecc, dens, dC.ecc); % interpolate dens. at radius ecc.
			d = dC; % use radius table
		case 'dend' % ganglion cell dendrite radius
			e = d.ecc; % functional eccentricity (deg)
			r = polyval(m.p.radGangCoef, e); % g.c. dend. radius (deg) *** use readTab
	end
	d.cover = pi * r .^ 2 .* dens; % coverage = area * density
	d.Properties.VariableDescriptions{'cover'} = 'Ganglion cell coverage';

function [d, m] = densCone(~, m) % cone density: Packer (89)

	% Load and adjust data table
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Packer (89) dens']; % data folder
	d = readFolder(folder); % create table from folder
	s = split(d.label, ", "); % split the label into quad, near
	d.quad = s(:, 1); d.near = s(:, 2); % store
	d.quad = categorical(d.quad); d.near = double(d.near); % fix the classes
	d = removevars(d, 'label'); % no longer required
	d = movevars(d, {'quad', 'near'},'before', 'ecc'); % match old order
	d = sortrows(d, {'quad', 'ecc'}); % sort rows
	
	% Convert and store
	s = 1 - .022; % linear shrinkage
	d.ecc = (1 / s) * d.ecc; % eccentricity, corrected for shrinkage (mm)
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	s = 1 - .05; % areal shrinkage
	dens = s * d.dens; % cone density, corrected for shrinkage (.001 x mm^-2)
	d.dens = 1000 * dens; % cone density (mm^-2)
	d.densDeg = (m.p.magRet ^ -2) * d.dens; % cone density (deg^-2)
	d.Properties.VariableDescriptions = {'Retinal quadrant', 'Near fovea', ...
		'Eccentricity (mm)', 'Cone density (mm^-^2)', 'Eccentricity (deg)', ...
		'Cone density (deg^-^2)'};

function [d, m] = densGang(~, m) % ganglion cell density: Wässle (89)

	% Load the empirical data, and convert mm to degrees
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Wässle (89)']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % rename label
	d.type = categorical(d.type); % for better listing
	d = sortrows(d, {'type', 'ecc'}); % sort by eccentricity
	d.dens = 10 .^ d.densLog; % density (mm^-2)
	eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	densDeg = (m.p.magRet) ^ -2 * d.dens; % density (deg^-2)

	% Correct for shrinkage, and store
	i = d.type == 'anatNarr'; % anatWide and func data are already corrected
	s = .9; % linear shrinkage
	eccDeg(i) = (1 / s) * eccDeg(i); % shift away from fovea
	densDeg(i) = s ^ 2 * densDeg(i); % reduce density
	d.eccDeg = eccDeg; d.densDeg = densDeg; % store
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'dens'} = 'Ganglion cell density (mm^-^2)';
	d.Properties.VariableDescriptions{'densDeg'} = ...
		'Ganglion cell density (deg^-^2)';

	%	Prepare for fitting
	d = sortrows(d, 'ecc'); % sort by eccentricity
	d.eccDegLog = log10(d.eccDeg); d.densDegLog = log10(d.densDeg); % for log fit 
	d.Properties.VariableDescriptions{'eccDegLog'} = 'Log eccentricity (deg)';
	d.Properties.VariableDescriptions{'densDegLog'} = ...
		'Log ganglion cell density (deg^-^2)';

function [d, m] = densGangFun(d, m) % g.c. dens., Wässle (89), vs. fun. ecc.

	% Initialise
	e = d.ecc; % functional eccentricity (deg): es x 1
	n = d.count; % ganglion cell count: es x 1
	
	% Differentiate cell count to find density versus functional eccentricity
	switch 'deriv' % method
		case 'deriv' % derivative
			switch 'all' % retinal region
				case 'near' % near fovea
					i = e < 10; % eccentricities close to fovea
					e = e(i); n = n(i); % select close eccentricities
			end
			es = 100; % number of points at which to resample e, for derivative
			eR = linspace(e(1), e(end), es); % resample eA at regular spacing
			n = interp1(e, n, eR); % interpolate to obtain density: es x 1
			dens = gradient(n, eR(2) - eR(1)) ./ (2 * pi * eR); % derivative
			eR(1) = 0; % extrapolate for first point
			dens(1) = interp1(eR(2: end), dens(2: end), 0, 'pchip'); % extrapolate
		case 'disc' % discrete
			dA = diff(pi * e .^ 2); % annulus areas (deg^2): (es - 1) x 1
			dn = diff(n); % count per annulus: (es - 1) x 1
			dens = (dn ./ dA); % density (deg^-2): (es - 1) x 1
			i = dn >= 2000; % annuli with at least 2000 cells
			e = e(i); es = length(e); % eccentricity (deg): es x 1
			dens = dens(i); % remove low-cell points
			dA = repmat(d(1, :), [es, 1]); % store
			dA.eccDeg = e; dA.densDeg = dens; % store new eccentricities, densities
	end
	
	% Store
	d = repmat(d(1, :), [es, 1]); % one row for each new eccentricity
	d.ecc = eR'; d.densDeg = dens'; % store new eccentricies, densities

function [d, m] = densGangSub(d, m) % split ganglion cells into subpopulations

	if isfield(m.densGangSub, 'z') % set name of variable to be split
		z = m.densGangSub.z; % user specified
	else, z = 'densDeg'; % default
	end
	switch 'fixed' % type of ratio
		case 'fixed' % ratio is constant across eccentricity
			d.type(:) = 'all'; % all ganglion cells
			dMid = d; % replicate
			dMid.type(:) = 'mid'; % midget ganglion cells
			r = m.p.ratGang; % ratio of midget to all ganglion cells
			dMid.(z) = r * dMid.(z); % reduce variable for midget ganglion cells
			d = [d; dMid]; % concatenate
		case 'var' % ratio varies with eccentricity
			dC = ratMidget(d, m); % ratio of midget to all ganglion cells
			r = interp1(dC.eccDeg, dC.ratio, d.eccDeg); % interpolate at ecc. in d
			
			% Calculate the densities of subpopulations
			dAll = d; dAll.type(:) = 'all'; % density of all ganglion cells
			dMid = d; dMid.type(:) = 'mid'; % midget ganglion cells
			dMid.densDeg = r .* dAll.densDeg; % density of midget g.c. (deg^-2)
			dOff = d; dOff.type(:) = 'off'; % off-centre midgets
			dOff.densDeg = m.p.ratOff * dMid.densDeg; % density of off-midgets (deg^-2)
			dOn = d; dOn.type(:) = 'on'; % on-centre midgets
			dOn.densDeg = (1 - m.p.ratOff) * dMid.densDeg; % den. of on-midgets (deg^-2)
			d = [dAll; dMid; dOff; dOn]; % concatenate
	end
	d.Properties.VariableDescriptions{'type'} = 'Cell type'; % describe

function [d, m] = doPrint(d, m) % print variable for storage in readTab

	v = d.(m.doPrint.var); % variable to print
	vs = length(v); % number of values to print
	n = 8; % number of values per line
	ls = ceil(vs / n); % number of lines to print
	for i = 1: ls % print a line
		j = n * (i - 1) + 1; % first index to print
		k = min(n * i, vs); % last index to print
		fprintf('%6.4f, ', v(j: k)); % print line
		fprintf('\n'); % finish line
	end

function [d, m] = eccFun(d, m) % calc. func. ecc.: McGregor (18), Schein (88)
		
	% Prepare McGregor (18) data
	[d1, m] = eccFunMcG(d, m); % load McGregor (18)
	if isfield(m.eccFunFit, 'mcGregor') % limit to reliable data
		i = d1.eccFun >= m.eccFunFit.mcGregor; % reliable functional eccentricities
		d1 = d1(i, :); % remove outlier data
	end

	% Prepare Schein (88) data
	[d2, m] = offset(d, m); % load Schein (88)
	i = d2.source == 'total'; d2 = d2(i, :); % use total offset (deg): es x 1
	eF = d2.ecc; % functional eccentricity (deg): es x 1
	o = d2.offset; % offset (deg): es x 1
	e = eF + o; % anatomical eccentricity (deg): es x 1
	d2.ecc = e; % anatomical eccentricity (deg): es x 1
	d2.eccFun = eF; % functional eccentricity (deg): es x 1
	if isfield(m.eccFunFit, 'schein') % limit Schein to high eccentricities
		i = d2.ecc >= m.eccFunFit.schein; % high eccentricities
		d2 = d2(i, :); % limit
	end

	% Combine
	d = outerjoin(d1, d2, 'mergeKeys', 1); % concatenate
	d = sortrows(d, {'source', 'ecc'}); % combine and sort

function [d, m] = eccFunFit(d, m) % fit func. ecc.: McGregor (18), Schein (88)

	% Fit functional eccentricity
	switch 'reg' % select method
		case 'hist' % histogram then spline

			% McGregor (18): bin eccentricity and average function ecc. in each bin
			i = d.source == 'mcGregor'; d1 = d(i, :); % select data
			bins = 10; % number of bins in which to average func. ecc.
			[~, edge, i] = histcounts(d1.ecc, bins); % find bin indices of func. ecc.
			binWid = edge(2) - edge(1); % bin width (deg)
			ecc = zeros(1, bins); eccFun = ecc; % allocate storage
			for j = 1: bins % loop over bins
				ecc(j) = edge(j) + .5 * binWid; % ecc. at middle of bin
				eccFun(j) = mean(d1.eccFun(i == j)); % average func. ecc. for this bin
			end
			d1 = d1(1: bins, :); % keep first few lines
			d1.ecc = ecc'; d1.eccFun = eccFun'; % store
		
			%	Interpolate all values with spline
			i = d.source == 'total'; % Schein values at high eccentricities
			d2 = d(i, :); % select data
			ecc = [d1.ecc; d2.ecc]; eccFun = [d1.eccFun; d2.eccFun]; % source data
			n = 50; % number of eccentricities at which to predict fun. ecc.
			eccInt = linspace(ecc(1), ecc(end), n)'; % interpolation ecc.
			eccFunInt = spline(ecc, eccFun, eccInt); % interpolate
			
		case 'reg' % regression

			% Fit polynomial to all data
			eccLog = log10(d.ecc); % for log versus log fit
			d.eccLog = eccLog; % store
			model = fitglm(d, 'poly3', 'predictorVars', 'eccLog', ...
				'responseVar', 'eccFun', 'link', 'log');

			% Rewrite data file
			n = 20; % number of rows
			source = repmat("model", [n, 1]); % data source
			eccLog = linspace(min(eccLog), max(eccLog), n)'; % log ecc. (deg)
			ecc = 10 .^ eccLog; % log anatomical eccentricity (deg)
			d = table(source, eccLog, ecc); % make table
			d.eccFun = predict(model, d); % functional eccentricity (deg)
			d.source = categorical(d.source); % for listing
			d.Properties.VariableDescriptions = {'Data source', ...
				'Log anatomical eccentricity (log deg)', ...
				'Anatomical eccentricity (deg)', 'Functional eccentricity (deg)'};
				% describe variables
			
			% Prepare to add points at high eccentricity
			eccInt = d.ecc; eccFunInt = d.eccFun; % match name with other methods
			dC = d(1: 2, :); d = [d; dC]; % add two rows to data file

		case 'spline' % regression then spline

			% Fit polynomial to McGregor (18)
			i = d.source == 'mcGregor'; d1 = d(i, :); % select data
			model = fitlm(d1, 'poly3', 'predictorVars', 'ecc', ...
				'responseVar', 'eccFun');
			eccEnd = d1.ecc(end); % save final eccentricity (deg)
			n = 20; % number of eccentricities at which to predict fun. ecc.
			d1 = d1(1: n, :); % keep first few lines
			eccBeg = 1; % assume x-intercept (deg)
			ecc = linspace(eccBeg, eccEnd, n)'; d1.ecc = ecc; % subsample ecc.
			eccFun = predict(model, d1); % calculate model predictions
			ecc(1) = eccBeg; eccFun(1) = 0; % intercept of prediction with ecc. axis
			d1 = d1(1: length(ecc), :); % trim to non-negative fun. ecc.
			d1.ecc = ecc; d1.eccFun = eccFun; % store
		
			%	Interpolate all values with spline
			i = d.source == 'total'; % Schein values at high eccentricities
			d2 = d(i, :); % select data
			ecc = [d1.ecc; d2.ecc]; eccFun = [d1.eccFun; d2.eccFun]; % source data
			n = 50; % number of eccentricities at which to predict fun. ecc.
			eccInt = linspace(eccBeg, ecc(end), n)'; % interpolation ecc.
			eccFunInt = spline(ecc, eccFun, eccInt); % interpolate

	end
	
	% Add points at high eccentricity
	u = m.p.magRet * 2.88; % eccentricity at which offset = 0 (deg):
		% from Schein (88) Figure 14 Nasal and Temporal
	eccInt = [eccInt; u; 90]; eccFunInt = [eccFunInt; u; 90];
		% add end-points (deg): es x 1
	d = d(1: n + 2, :); d.ecc = eccInt; d.eccFun = eccFunInt; % store
	d.eccLog = log10(d.ecc); % for log versus log fit

function [d, m] = eccFunMcG(~, m) % functional vs. anat. ecc.: McGregor (18)

	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'McGregor (18)']; % data folder
	d = readFolder(folder); % create table from folder
	d.eccMm = d.eccUm / 1000; d.eccFunMm = d.eccFunUm / 1000; % convert to mm
	d = renamevars(d, 'label', 'source'); % rename label
	d = sortrows(d, {'source', 'eccMm'}); % sort by source, eccentricity
	d.source = categorical(d.source); % make it categorical
	d.ecc = m.p.magRet * d.eccMm; d.eccFun = m.p.magRet * d.eccFunMm; % mm to deg
	d.Properties.VariableDescriptions{'source'} = 'Source';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';

function [d, m] = freq(d, m) % calculate spat. freq. cutoff versus eccentricity

	%	Calculate radius of excitatory cell centre mechanism
	import prim.readTab % find function
	e = [0, 1, 3, 10, 30]; % eccentricities commonly used in runPrim (deg): 1 x es
	rOpt = polyval(m.p.radOptCoef, e); % radius of point spread function (deg)
	rGang = readTab('radGang', e); % radius of ganglion cell dendritic tree (deg)
	fac = readTab('ratCort', e); % inverse cort. mag'n factor (deg/mm)
	rBeta = fac * m.p.radBetaMm; % 4CBeta convergence radius (deg)
	rCen = sqrt(rOpt .^ 2 + rGang .^ 2 + rBeta .^ 2); % ex. cell radius (deg)

	%	Calculate and display spatial frequency cutoff for excitatory cell
	a = m.freq.cutoff; % centre mechanism attenuation at cutoff
	freq = sqrt(- log(a)) ./ (pi * rCen); % spatial frequency cutoff (cycles/deg)
	fprintf('%6.4f, ', freq); fprintf('\n'); % print separations

function [d, m] = hist(d, m) % calculate histogram: data can be in a structure

	% Initialise
	if isfield(m.hist, 'edges') % set edges or number of bins
		arg = m.hist.edges; % user-specified
	else
		arg = 10; % number of bins
	end
	
	%	Obtain x data
	x = string(m.x); % location of x data: scalar string
	x = split(x, "."); % split m.x at dots: string array
	xs = length(x); % number of substrings
	if xs > 1 % data are in a structure
		x = d.(x(1)).(x(2)); % x data: should be numeric
	else % there is only one string
		x = d.(x); % x data
	end

	% Calculate histogram
	[n, e] = histcounts(x, arg); % obtain histogram counts and edges
	c = e(1: end - 1) + .5 * diff(e); % turn edges into bin centres: 1 x bs

	%	Store
	d = d(1, :); % keep only the first row
	d.cen = c; % bin centres
	d.edge = e; % bin edges
	d.(m.y) = n; % counts
	d.Properties.VariableDescriptions{'cen'} = 'Bin centre';
	d.Properties.VariableDescriptions{'edge'} = 'Bin edge';
	d.Properties.VariableDescriptions{m.y} = 'Count';

function [d, m] = mtf(~, m) % compile modulation transfer function

	% Parameters for modulation transfer function: Navarro (93)
	ecc = [	0			5			10		20		30		40		50		60]'; % eccentricity (deg)
	a =		[	.172	.245	.245	.328	.606	.82		.93		1.89]';
	b =		[	.037	.041	.041	.038	.064	.064	.059	.108]';
	c =		[	.22		.2		.2		.14		.12		.09		.067	.05]';
	
	% Calculate MTF and store
	es = length(ecc); % number of eccentricity
	freq = linspace(0, 60); % spatial frequency (cycles/deg): 1 x fs
	mtf = (1 - c) .* exp(- a .* freq) + c .* exp(- b .* freq); % MTF: es x fs
	freq = repmat(freq, [es, 1]); % spatial frequency (cycles/deg): es x fs
	d = table(ecc, a, b, c, freq, mtf); % store as table: es x 6
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'freq'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'mtf'} = 'Modulation transfer function';

function [d, m] = offset(~, m) % ganglion cell lateral offset: Schein (88)

	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Schein (88)']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d = sortrows(d, {'source', 'eccMm'}); % sort by source, eccentricity
	d.source = categorical(d.source); % make it categorical
	d.ecc = m.p.magRet * d.eccMm; d.offset = m.p.magRet * d.offsetMm; % mm to deg
	d.Properties.VariableDescriptions{'source'} = 'Source';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'offset'} = 'Displacement (deg)';

function [d, m] = prep(d, m) % prepare for plot

	% Set independent and dependent variables, and list preparations
	if isfield(m.prep, 'x'), nX = m.prep.x; x = d.(nX); % independent variable
	elseif isfield(m, 'x'), nX = m.x; x = d.(nX); end
	%{
	if isfield(m.prep, 'y'), nY = m.prep.y; y = d.(nY); % dependent variable
	elseif isfield(m, 'y'), nY = m.y; y = d.(nY); end
	%}
	n = fieldnames(m.prep)'; % names, including preparations (cell): 1 x ps

	%	Prepare
	for p = string(n) % loop over preparations
		switch p % choose preparation
			case 'linspace' % replace x with linearly spaced vector *** keep?
				x = linspace(m.prep.linspace{:}); % vector
				d.(nX) = x(:); % store as column vector
			case 'polar' % turn real and imaginary parts into amp., phase *** delete
				i = d.real == 1; % rows containing real parts
				dC = d(~ i, :); d = d(i, :); % tables with real, imaginary parts
				resp = complex(d.resp, dC.resp); % combine to obtain complex response
				d.amp = abs(resp); % amplitude (Hz)
				phase = angle(resp); % phase (radians)
				phase = unwrap(phase); % avoid 360 deg jumps in phase
				d.phase = (180 / pi) * phase; % phase (deg)
			case 'wool' % components to amplitude, phase
		
				% Convert components to amplitude, phase
				x = shiftdim(d.(m.x), 1); % predictor variables: rs x vs
				i = logical(x(:, 1)); % component, 1 for cosine: rs x 1
				c = x(i, [2, 3]); % contrast: os x 2
				f = x(i, 4); % spatial frequency (cycles / deg): os x 1
				r = d.(m.y); % predicted response: 1 x rs
				r = r(i) + 1i * r(~ i); % change cosine and sine comp. to complex: 1 x os
				act = abs(r); % amplitude (Hz): 1 x os
				phase = 180 * angle(r) / pi; % phase (deg): 1 x os
		
				% Store
				os = size(c, 1); % number of observations
				d = repmat(d, [os, 1]); % one row per observation: os x 1
				d.cont = c; % contrast: conts x 2
				d.freq = f; % spatial frequency (cycles / deg): os x 1
				d.act = act'; % amplitude (Hz): os x 1
				d.phase = phase'; % phase (deg): 1 x os
				m.x = 'freq'; % update the independent variable
	
			case 'zero' % replace 0 on log x-axis by positive value *** check ***
				i = x == 0; % rows for which x is 0
				x(i) = m.prep.zero; % replace
				d.(nX) = x; % store
		end
	end

function [d, m] = priv(d, m) % enable private line from parafoveal cone to g.c

	import prim.readTab; % function location
	switch m.priv.var % densGang or radGang
		case 'densGang' % set parafoveal off- and on-g.c dens. equal to cone dens.

			%	Initialise
			e = [0: .25: .75, 1: .5: 2.5, 3: 1: 9, 10: 5: 50]; % ecc. (deg): 1 x es
			densCone = readTab('densCone', e); % cone density (cells/deg^2): 1 x es
			dens = interp1(d.ecc, d.densDeg, e); % g.c. dens. (cells/deg^2): 1 x gs
			i = find(e == 1); j = find(e == 3); % transition points
			k = [1: i, j: length(e)]; % empirical points
			l = i + 1: j - 1; % points to interpolate

			% Construct densities and store
			densOff = densCone; % off-centre ganglion cell density = cone density
			densOff(j: end) = m.p.ratSign * dens(j: end); % off-cell density
			densOff(l) = spline(e(k), densOff(k), e(l)); % interpolate
			densOn = densCone; % on-centre ganglion cell density = cone density
			densOn(j: end) = (1 - m.p.ratSign) * dens(j: end); % on-cell density
			densOn(l) = spline(e(k), densOn(k), e(l)); % interpolate
			d = table(densOff, densOn); % store
			
		case 'radGang' % reduce parafoveal ganglion cell dendritic radius

			%	Generate reduced vector of eccentricity
			e = [0, 1, 9: 3: d.eccFun(end)]; % ecc. in steps of 1, 3 deg: 1 x es
			es = length(e); % number of eccentricities
			d = repmat(d(1, :), [es, 1]); % reduce table to match vector length
			d.eccFun = e'; % new eccentricities (deg): 1 x es
		
			%	Calculate ganglion cell dendritic radii, and store
			r = predict(m.model{1}, d); % new radii (deg): 1 x es
			r(1: 2) = 0; % private line defined by zero radius
			%	r(2: 6) = linspace(r(2), r(6), 5); % interpolate between 0 and empirical
			d.radius = r; % store in d
			%	m.p.radGangTab = [d.eccFun, d.radius]'; % for storage in readTab

	end

function [d, m] = procIm(~, m) % process ImageNet images: Deng (09)

	% Initialise
	import prim.filtIm prim.readIm prim.readImFolder; % import functions
	if isfield(m.procIm, 'files') % serial numbers of image files
		file = m.procIm.files; % user-specified
	else
		file = []; % default is all
	end
	proc = string(m.procIm.proc); % screen or filter
	widMax = 500; % maximum width and height of stored image (pix): 1 x 1

	%	Initialise table
	[folder, name] = readImFolder(file); % list names of image files
	fs = length(name); % number of files
	if isempty(file), file = 1: fs; end % serial numbers of all image files
	d = table(nan, "", nan(1, 3), false, nan, nan, ...
		'variableNames', {'file', 'name', 'sizeOrig', 'valid', 'width', 'waveDom'});
	d.wave = struct('wave', nan); % preferred wavelength: structure
	d.image = struct('im', nan); % processed image: structure
	d = repmat(d, [fs, 1]); % one row for each file
	d.Properties.VariableDescriptions = ...
		{'Image file number', 'File name', 'Original image size (pixels)', ...
		'Image validity', 'Image width (pixels)', ...
		'Dominant wavelength (pixels)', 'Preferred wavelength (pixels)', ...
		'Image (sRGB-units)'};

	% Read, process, and store images
	%	for i = 1: fs % loop over files
	parfor i = 1: fs % loop over files using parallel processing

		%	Read file
		dC = d(i, :); % current row of table
		nameC = name(i); % name of current file
		[~, nameN, ~] = fileparts(nameC); % name with no extension
		dC.file = file(i); dC.name = nameN; % store
		[im, s] = readIm(folder, nameC, widMax); % read and crop image file
		if ~ isempty(im) % file is readable and three-dimensional

			%	This is a valid image
			dC.sizeOrig = s; % original size of image (pix): [ys, xs, cs]
			wid = size(im, 1); % width of current image (pix)
			dC.valid = true; dC.width = wid; % store
			dC.image = struct('im', im); % image cropped to acceptable size and ...
				% square: ys x xs x cs, uint8
			
			if proc == "filter" % screening is finished, now filter image

				%	Filter image to find dominant wavelength
				[mag, t] = filtIm(im); % Gabor magnitude: % ys x xs x gs, gs = ws * dirs
				[~, j] = max(mag, [], 3); % Gabor with maximum magnitude: ys x xs
				w = t(j(:), 1); % preferred wavelength (pix): ys * xs x 1
				waveP = uint16(reshape(w, size(j)));
					% preferred wavelength (pix): ys x xs, uint16 for reduced storage size
				dC.wave = struct('wave', waveP); % store preferred wavelength
	
				%	Calculate dominant wavelength
				%	wMax = quantile(w, .9); % wavelength at 90% decile (pix):
				w = w(w < max(w)); % wavelengths less than maximum
					% wavelengths in top bin are typically from the background
				dC.waveDom = median(w); % store dominant wavelength (pix)
				if isnan(dC.waveDom) % this occurs when all pref. wavelengths are equal
					dC.valid = false; % declare results invalid
				end

			end

		end
		d(i, :) = dC; % store current row in table
		
	end

function [d, m] = psf(d, m) % calculate point spread function

	% Calculate point spread function as inverse Fourier transform of MTF
	es = size(d, 1); % number of eccentricities
	a = d.a; b = d.b; c = d.c; % coefficients: es x 1
	w = .2; xs = 100; % location range (deg), number of locations
	x = linspace(-.5, .5, xs + 1) * w; % location (deg): 1 x (xs + 1)
	x(end) = []; % make it open-ended: 1 x xs
	y = 2 * a .* (1 - c) ./ (a .^ 2 + (2 * pi * x) .^ 2) + ...
		2 * b .* c ./ (b .^ 2 + (2 * pi * x) .^ 2); % PSF: es x xs

	% Fourier transform PSF to check on solution
	f = linspace(- .5, .5, xs + 1) * (xs / w); % spatial freq. (cycles/deg)
	f(end) = []; % make it open-ended: 1 x xs
	f = fftshift(f); % shift negative frequencies to end: 1 x xs
	z = fft(y, [], 2) ./ sum(y, 2); % Fourier transform: es x xs
		% why are odd components negative?
	z = abs(z); % MTF: es x xs

	% Store in table
	d.loc = repmat(x, [es, 1]); % location (deg): es x xs
	d.psf = y; % PSF: es x xs
	d.Properties.VariableDescriptions{'loc'} = 'Location (deg)'; % describe
	d.Properties.VariableDescriptions{'psf'} = 'Point spread function';
	d.freq = repmat(f, [es, 1]); % spatial frequency (cycles/deg): es x xs
	d.mtf = z; % MTF: es x xs
	d.Properties.VariableDescriptions{'freq'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'mtf'} = 'Modulation transfer function';

function [d, m] = rad(d, m) % calculate radii contibuting to centre, surround

	% Calculate optical point spread function radius
	d = radOpt(d, m); % load optical data
	switch 'fun' % choose method
		case 'emp' % use eccentricities in optics empirical data
			d = d(d.ecc >= 0, :); % retain only temporal data
			eF = d.ecc; % eccentricity (deg)
			d.optics = d.radius; % PSF radius (deg)
		case 'fun' % calculate eccentricities for fitted data
			l = log10([.006, 40]); % log eccentricity limits (deg)
			es = 30; % number of eccentricities
			eF = logspace(l(1), l(2), es)'; % eccentricity (deg): es x 1
			r = polyval(m.p.radOptCoef, eF);
			d = repmat(d(1, :), [es, 1]); % repeat first row of table
			d.ecc = eF; d.optics = r; % store
	end

	% Calculate remaining radii
	import prim.readTab; % function location
	d.fieldHor = exp(polyval(m.p.radHorCoef, eF)); % receptive field radius (deg)
	%	d.denGang = polyval(m.p.radGangCoef, eF); % ganglion cell dend. radius (deg)
	d.denGang = readTab('radGang', eF, m); % ganglion cell dend. radius (deg)
	d.cenConv = sqrt(d.optics .^ 2 + d.denGang .^ 2); % centre radius (deg)
	d.surConv = sqrt(d.optics .^ 2 + 2 * d.fieldHor .^ 2); % sur. radius (deg)

	% Turn lateral array into vertical
	s = {'optics', 'fieldHor', 'denGang', 'cenConv', 'surConv'};
		% sources of radius (deg)
	ss = length(s); % number of sources
	dC = cell(ss, 1); % allocate storage
	for i = 1: ss % loop over sources
		dCC = d; % replicate existing table
		sC = s{i}; % current source
		dCC.source(:) = sC; dCC.radius = d.(sC); % current source and radius
		dC{i} = dCC; % store
	end
	d = vertcat(dC{:}); % concatenate into vertical array
	d = removevars(d, ['rrfLog', 'rrf', s]); % variables not needed
	d.source = categorical(d.source); % make it categorical, for listing purposes
	d.radiusMin = 60 * d.radius; % convert to minutes, for plotting
	d.Properties.VariableDescriptions = {'Source of radius', ...
		'Eccentricity (deg)', 'Radius (deg)', 'Radius (min)'};

function [d, m] = radBeta(~, m) % 4CBeta convergence radius: multiple sources

%	Sources: Blasdel (83), Freund (89), Yabuta (98), Bauer (99): details in file d

	% Load data
	load([m.folder, filesep, '4CBeta radius/gencortrads.mat'], 'd');
	%	d = gencortrads; % table
	d = renamevars(d, {'Radii', 'Location', 'Source'}, ...
		{'radius', 'loc', 'source'}); % rename for convenience
	d.loc = categorical(d.loc); d.source = categorical(d.source); % for plotting
	d.Properties.VariableDescriptions = {'Radius (mm)', 'Stage', 'Source'};

	% Calculate convergence radius
	i = d.loc == 'LGN axon terminal'; % indices
	radGen = mean(d.radius(i)); % mean geniculate radius (mm)
	%	radBeta = mean(d.radius(~ i)); % mean 4CBeta dendritic radius (mm)
	i = d.source == 'Bauer et al., 1999'; % refers to 25 cells in Lund (80)
	radBeta = d.radius(i); % 4CBeta dendritic radius (mm)
	radConv = sqrt(radGen ^ 2 + radBeta ^ 2); % cross-correlation (mm)
	d(end + 1, :) = d(1, :); % add a row
	d.radius(end) = radConv; d.loc(end) = 'Convergence'; d.source(end) = 'All';

function [d, m] = radBetaEcc(d, m) % geniculocortical convergence radius

	d.radius = d.fac * m.p.radBetaMm; % radius (deg) as a function of eccentricity
	d.Properties.VariableDescriptions{'radius'} = 'Convergence radius (deg)';

function [d, m] = radCen(~, m) % compile ganglion cell centre radius data

	% Croner (95)
	load([m.folder, filesep, 'Croner (95) rad/parvo_cenrad_ecc.mat'], 'b');
	ecc = b(:, 1); % eccentricity (deg): es x 1
	es = size(ecc, 1); % number of eccentricities
	radius = b(:, 2); % ganglion cell centre radius (deg): es x 1
	source = repmat("Croner", [es, 1]); % data source
	d1 = table(source, ecc, radius); % store

	% Lee (98)
	load([m.folder, filesep, 'Lee (98)/CenLeeexp.mat'], 'CenLeeexp');
	d2 = CenLeeexp; % ecc. and radius: es x 2
	ecc = d2(:, 1); % eccentricity (deg): es x 1
	es = size(ecc, 1); % number of eccentricities
	dev = 10 .^ d2(:, 2); % standard deviation (min): es x 1
	radius = sqrt(2) * dev / 60; % gang. cell centre radius (deg): es x 1
	source = repmat("Lee", [es, 1]); % data source
	d2 = table(source, ecc, radius); % store
	
	%	Godat (22)
	d3 = radGangField(d2, m); % read data
	i = d3.type == 'centre'; d3 = d3(i, :); % keep only centre data
	d3 = renamevars(d3, 'type', 'source'); % prepare to set source
	d3.source(:) = "Godat"; d3.source = string(d3.source); % set source
	d3 = removevars(d3, {'eccLog', 'radLog'}); % remove extraneous variables

	% Concatenate
	d = [d1; d2; d3]; % concatenate
	d.source = categorical(d.source); % for listing
	d.radiusLog = log10(d.radius); % add logarithm for fitting purposes
	d.radiusMin = 60 * d.radius; % radius (min) for display purposes
	d.dev(:) = nan; d.dev(end - es + 1: end) = dev; % s.d. (min)
	d.Properties.VariableDescriptions = {'Source', ...
		'Eccentricity (deg)', 'Centre radius (deg)', ...
		'Log of centre radius (deg)', 'Centre radius (min)', ...
		'Standard deviation (min)'};

function [d, m] = radCone(~, m) % cone inner segment radius: Packer (89) diam

	% Generate a table for each quadrant
	folder = [m.folder, filesep, 'Packer (89) rad']; % data folder
	f = dir([folder, '/*.mat']); % mat-file list: fs x 1
	fs = length(f); % number of files
	dC = cell(fs, 1); % allocate storage			
	for i = 1: fs % loop over files

		% Collect data
		fC = f(i).name; % full name of current file
		q = fC(1: 2); % name of current quadrant: 2 x 1
		quad = categorical(string(q)); % store
		d = load([folder, '/', fC]); % contents of file
		d = d.(q); % eccentricity and diameter (mm, um): ds x 2
		ds = size(d, 1); % number of rows

		% Store in table
		quad = repmat(quad, [ds, 1]); % match number of rows
		dC{i} = table(quad, d(:, 1), d(:, 2), 'variableNames', ...
			{'quad', 'ecc (mm)', 'diam (um)'});

	end

	% Combine and store
	d = vertcat(dC{:}); % combine data tables
	d = sortrows(d, {'quad', 'ecc (mm)'}); % sort rows
	d.Properties.VariableDescriptions = ...
		{'Retinal quadrant', 'Eccentricity (mm)', 'Cone diameter (um)'};

function [d, m] = radHor(d, m) % hor. cell radius: Wässle (89), Packer (02)

	% Dendrite radius: Wässle (89) Horizontal
	[d, m] = radHorDend(d, m); % obtain data
	d = renamevars(d, 'type', 'source'); % rename label
	i = d.source == "H1"; d = d(i, :); % keep only H1 cells
	d.source(:) = "dend"; % source of data
	d = d(d.eccFun <= 16, :); % field only at low eccentricities
	r = d.area; % dendrite area (mm^2)
	r = sqrt(r / pi); % radius (mm): es x 1
	d.radius = m.p.magRet * r; % radius (deg): es x 1
	d.Properties.VariableDescriptions{'radius'} = 'Radius (deg)';
	d1 = d; % store
	
	% Field radius: Packer (02)
	[d, m] = radHorField(d, m); % obtain data
	d = outerjoin(d1, d, 'mergeKeys', 1); % concatenate

function [d, m] = radHorDend(~, m) % hor. cell dendrite area: Wässle (89) Hor

	% Compile data into table
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Wässle (89) Horizontal']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % rename label
	d.type = categorical(d.type); % make it categorical
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	d.area = 10 .^ d.areaLog; % store the anti-log
	d.Properties.VariableDescriptions{'type'} = 'Cell type'; % describe variable
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'areaLog'} = 'Log area (mm^2)';
	d.Properties.VariableDescriptions{'area'} = 'Area (mm^2)';

	% Add functional eccentricity
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	d.eccFun = interp1(dE.ecc, dE.eccFun, d.eccDeg); % interpolate on fun. ecc.
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';

function [d, m] = radHorField(~, m) % horizontal cell field radius: Packer (02)

	% Compile data into table
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Packer (02)']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d.source = categorical(d.source); % make it categorical
	d.eccFun = m.p.magRet * d.ecc; % functional eccentricity: (deg)
	r = d.diam; % diameter (um)
	r = .5 * m.p.magRet * .001 * r; % radius (deg)
	d.radius = r / sqrt(log(10)); % convert from .1 * peak to 1 / e (deg)
	d.Properties.VariableDescriptions{'source'} = 'Source of radius';
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Field radius (deg)';
	
	% Bin the data, split bins into quarters, convolve radius of each quarter
	i = d.source == 'wide' & d.eccFun >= 16 & d.eccFun < 60; % relevant rows
	dC = d(i, :); % select wide fields and eccentricities with sufficient data
	bins = 8; % number of bins
	binRads = 4; % number of radius bins
	[~, edge, bin] = histcounts(dC.eccFun, bins); % ecc. bin edges and indices
	rConv = zeros(1, binRads); % allocate storage
	for i = 1: bins % loop over bins
		dCC = dC(bin == i, :); % rows for this bin
		r = dCC.radius; % radii in this bin (deg)
		[~, ~, binRad] = histcounts(r, binRads); % bin number for radii
		rC = zeros(1, binRads); % allocate storage
		for j = 1: binRads % loop over radius bins
			rC(j) = mean(r(binRad == j)); % mean of radii in current radius bin (deg)			
		end
		rConv(i) = sqrt(sum(rC .^ 2, 'omitnan')); % convolution of radii (deg)
	end
	
	% Store
	dC = repmat(dC(1, :), [bins, 1]); % one row for each bin
	eccCen = edge(1: bins) + .5 * (edge(2) - edge(1)); % centres of ecc. bins
	dC.eccFun = eccCen'; % ecc. bin centres (deg)
	dC.source(:) = 'conv'; % source is convoluted radius
	dC.radius = rConv'; % convoluted radius (deg)
	d = [d; dC]; % concatenate

function [d, m] = radGang(~, m) % ganglion cell denritic radius: Watanabe (89)

	% Load radius data and convert eccentricity
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Watanabe (89)']; % data folder
	d = readFolder(folder); % create table from folder
	d.label = []; % remove this variable
	d.ecc = m.p.magRet * d.eccMm; % eccentricity (deg): es x 1
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	d.eccFun = interp1(dE.ecc, dE.eccFun, d.ecc); % convert to functional ecc.
	d = sortrows(d, 'eccFun'); % sort on eccentricity
	%	e = d.eccFun; % functional eccentricity (deg)
	diam = 10 .^ d.diamLog; % dendritic field diameter (um): es x 1
	r = (.5 * diam / 1000) * m.p.magRet; % radius (deg): es x 1
	
	%	Reduce radius at low eccentricities to ensure a cone-g.c. private line
	%{
	i = 1: 8; is = length(i); % eccentricities below 3 deg
	e1 = 0; r1 = .001; % radius at lowest eccentricity
	r(i) = r1 + (r(is) - r1) * (e(i) - e1) / (e(is) - e1); % interpolate
	%}

	% Store
	d.diam = diam; d.radius = r; % store
	d.Properties.VariableDescriptions = ...
		{'Eccentricity (mm)', 'Log of dendritic field diameter (um)', ...
			'Eccentricity (deg)', 'Functional eccentricity (deg)', ...
			'Dendritic field diameter (um)', 'Dendritic field radius (deg)'};

function [d, m] = radGangField(~, m) % g.c. centre, surround radius: Godat (22)

	% Load data and create a table
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Godat (22)']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % remove this variable
	d.type = categorical(d.type); % make it categorical
	d.ecc = 10 .^ d.eccLog; d.radius = 10 .^ d.radLog; % take anti-logs
	d = sortrows(d, {'type', 'ecc'}); % sort by type, eccentricity
	d.Properties.VariableDescriptions{'type'} = 'Mechanism';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Mechanism radius (deg)';

function [d, m] = radOpt(~, m) % point spread func.: Navarro (93), Williams (81)

	% Generate a table from the data
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Navarro (93)']; % data folder
	d = readFolder(folder); % create table from folder
	d = sortrows(d, {'label', 'ecc'}); % sort on version, eccentricity
	d = renamevars(d, 'label', 'source'); % rename label
	d.source(:) = "optics"; d.source = categorical(d.source); % clean label

	% Clean data
	r = 10 .^ d.rrfLog; % retinal resolution function (min)
	es = .5 * size(d, 1); % number of eccentricities
	r = mean([r(1: es), r(es + (1: es))], 2); % mean over versions
	d = d(1: es, :); % keep only mean
	d.ecc = [-60, -50, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 50, 60]';
		% clean eccentricities (deg)
	d.rrf = r; % retinal resolution function, given as	full width
		%	at half-height of Gaussian (min): es x 1

	% Convert rrf into radius
	r = r / (2 * sqrt(log(2))); % radius of point spread function (min): es x 1
	r = r / 60; % radius of point spread function (deg): es x 1
	d.radius = m.p.ratOpt * r; % convert human to macaque (deg)
		%	Williams (81) Figure 5a, 13 weeks old, 6 mm pupil, LSF: 1/e radius = 1 min
	d.Properties.VariableDescriptions{'source'} = 'Source of data';
	d.Properties.VariableDescriptions{'rrf'} = 'Retinal resolution function (min)';
	d.Properties.VariableDescriptions{'radius'} = 'PSF radius (deg)';

function [d, m] = radSur(~, m) % g.c. sur.: Croner (95), Lee (98), Godat (22)

	% Compile data into table: Croner (95)
	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Croner (95) rad']; % data folder
	d = readFolder(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d.source = categorical(d.source); % make it categorical
	d.ecc = 10 .^ d.eccLog; % eccentricity (deg)
	d.radius = 10 .^ d.radiusLog; % surround radius (deg)
	d.Properties.VariableDescriptions{'source'} = 'Source of radius'; % describe
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Surround radius (deg)';

	% Lee (98): load both centre and surround data to find surround eccentricities
	folder = [m.folder, filesep, 'Lee (98)']; % data folder
	dC = readFolder(folder); % centre s.d. data
	dC = renamevars(dC, 'label', 'source'); % rename label
	dC.source = categorical(dC.source); % make it categorical
	dC.radCen = 10 .^ dC.radLog; % centre s.d. (min)
	dC = sortrows(dC, 'radCen'); % sort to match with surround data
	load([folder, filesep, 'CenSurLee'], 'CenSurLee'); % load surround
	dS = CenSurLee; % surround data
	dS = table(dS(:, 1), dS(:, 2), 'variableNames', {'radCen', 'radSur'}); % table
	dS = sortrows(dS, 'radCen'); % sort to match with centre data
	dC.radius = sqrt(2) * dS.radSur / 60; % radius (deg)
	dC.source(:) = "Lee"; % show source
	d = outerjoin(d, dC, 'mergeKeys', 1); % combine Croner, Lee
	
	% Godat (22)
	dC = radGangField(d, m); % read data
	i = dC.type == 'surround'; dC = dC(i, :); % keep only surround data
	dC = renamevars(dC, 'type', 'source'); % prepare to set source
	dC.source(:) = "Godat"; dC.source = categorical(dC.source); % set source
	d = outerjoin(d, dC, 'mergeKeys', 1); % add Godat

	%	Combine
	d.radiusMin = 60 * d.radius; % convert to minutes, for plotting

function [d, m] = ratCort(~, m) % inverse cort. magnification factor: Dow (81)

	% Calculate factor
	if isfield(m.ratCort, 'eccs') % number of eccentricities
		es = m.ratCort.eccs; % user-defined
	else % 
		es = 20; % default
	end
	%{
	xLim = log10([7, 3000]); % limits for log x 
	eMin = logspace(xLim(1), xLim(2), es)'; % eccentricity (min): es x 1
	%}
	eMin = linspace(0, 3000, es)'; % eccentricity (min): es x 1
	fMin = log10(eMin) - 1.5; % mysterious conversion
	fMin = .8124 + .5324 * fMin + .0648 * fMin .^ 2 + .0788 * fMin .^ 3; % fit
	fMin = 10 .^ fMin; % inverse CMF (min/mm): es x 1
	e = eMin / 60; % eccentricity (deg): es x 1
	eLog = log10(e); % log eccentricity (deg), for fitting: es x 1
	f = fMin / 60; % inverse CMF (deg/mm): es x 1
	w = (f(1) ./ f) .^ 2; % regression weights, high for low eccenticity: es x 1
	
	% Store
	s = repmat("dow", [es, 1]); % data source
	d = table(s, e, eMin, eLog, f, fMin, w, 'variableNames', ...
		{'source', 'ecc', 'eccMin', 'eccLog', 'fac', 'facMin', 'weight'}); % table
	d.Properties.VariableDescriptions = ...
		{'Data source', 'Eccentricity (deg)', 'Eccentricity (min)', ...
		'Log eccentricity (deg)', 'ICMF (deg/mm)', 'ICMF (min/mm)', ...
		'Regression weight'}; % add descriptions
	d.source = categorical(d.source); % for listing purposes

	% Add foveal representation
	%{
	d0 = d(1, :); % first row only
	d0.source = "zero"; % distinguish from published data
	d0.ecc = 0; % zero eccentricity
	d0.fac = interp1(d.ecc, d.fac, 0, 'spline');
	d = [d; d0]; % contenate data tables
	d = sortrows(d, 'ecc'); % sort by eccentricity
	%}
	fac = interp1(d.ecc(2: end), d.fac(2: end), 0, 'spline'); % extrapolate
	d.ecc(1) = 0; d.fac(1) = fac; % store

function [d, m] = ratCount(d, m) % ratio of cone to ganglion cell count

	eInc = .5; % eccentricity increment (deg)
	e = (0: eInc: d.ecc(end))'; % eccentricities at which to interpolate
	cG = interp1(d.ecc, d.count, e); % ganglion cell count
	dC = densCone(d, m); dC = countCone(dC, m); % cone count
	cC = interp1(dC.eccDeg, dC.count, e); % ganglion cell count
	r = cC ./ cG; % ratio of cone to ganglion cell count
	d = repmat(d(1, :), [length(e), 1]); % resize table
	d.ecc = e; d.ratio = r; % store eccentricity, ratio
	d.densDeg = cG ./ (pi * e .^ 2); % cumulative ganglion cell density (deg^-2)
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Ratio of cone to ganglion cell count'; % describe

function [d, m] = ratDens(d, m) % ratio of cone to ganglion cell density

	dC = densCone(d, m); % load cone density
	i = dC.quad == 'temporal'; dC = dC(i, :); % keep only temporal quadrant
	densC = interp1(dC.eccDeg, dC.densDeg, d.ecc); % cone dens. at ecc. in d
	d.ratio = densC ./ d.densDeg; % ratio of cone to ganglion cell density
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Ratio of cone to ganglion cell density'; % describe
	
function [d, m] = ratMidget(~, m) % mid. / all g.c.: % Dacey (94), Grünert (93)

	import prim.readFolder; % function is located in namespace prim
	folder = [m.folder, filesep, 'Dacey (94)']; % data folder
	d = readFolder(folder); % log ratio versus eccentricity
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	r = 10 .^ d.ratioLog; % ratio of midget ganglion cells to all g.c.
	d.ratioPer = min(100 - 6.5, r); ...
		% Grünert et al. found that 5-8% of fovel ganglion cells were parasol
	d.ratio = .01 * d.ratioPer; % midget/ all g.c.
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'ratioPer'} = ...
		'Midget / all ganglion cells (%)';
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Midget / all ganglion cells';

function [d, m] = resp(~, m) % g.c. spatial frequency resp.: multiple authors

	import prim.readFolder; % function is located in namespace prim
	switch m.resp.source % data source
		case 'Croner' % Croner (95)
			folder = [m.folder, filesep, 'Croner (95) resp'];
			d = readFolder(folder); % impulse rate versus spatial frequency
			d = renamevars(d, 'label', 'source'); % rename for consistency
			d.source = categorical(d.source); % for pretty printing
			d.freqS = 10 .^ d.freqLog; % spatial frequency (cycles/deg)
			d.freqT(:) = 4.22; % temporal frequency (Hz)
			d.amp = 10 .^ d.respLog; % response (Hz)
			d.phase(:) = 0; d.resp = d.amp; % store
			d.Properties.VariableDescriptions{'amp'} = 'Response amplitude (Hz)';
		case 'Wool' % Wool (18)
			f = [m.folder, filesep, 'Wool (18)/']; % data folder
			load([f, 'R1.mat'], 'R1'); % log amplitude (Hz): fs x 2
			load([f, 'Rphase.mat'], 'Rphase'); % phase (deg): fs x 2
			fs = size(R1, 1); % number of frequencies
			source = repmat("Wool", [fs, 1]); source = categorical(source); % source
			d = table(source); % make output table
			R1 = sortrows(R1); Rphase = sortrows(Rphase); % sort on frequency: fs x 2
			f = mean([R1(:, 1), Rphase(:, 1)], 2); % log frequency (cyc/deg): fs x 1
			d.freqS = 10 .^ f; % spatial frequency (cycles/deg): fs x 1
			d.freqT(:) = 2; % temporal frequency (Hz)
			a = 10 .^ R1(:, 2); % response amplitude (Hz): fs x 1
			p = Rphase(:, 2); % response phase (deg): fs x 1
			r = a .* exp(1i * pi * p / 180); % response (Hz), complex: fs x 1
			d.amp = a; d.phase = p; d.resp = r; % store
			d.Properties.VariableDescriptions{'amp'} = 'Response amplitude (Hz)';
			d.Properties.VariableDescriptions{'phase'} = 'Response phase (deg)';
		case 'Yeh' % Yeh (95)

			% Sensitivity data
			folder = [m.folder, filesep, 'Yeh (95) sens'];
			d = m.readFolder(folder); % contrast sensitivity versus temporal frequency
			d = renamevars(d, 'label', 'source'); % rename for consistency
			d.source = categorical(d.source); % for listing
			d.freqS(:) = 0; % spatial frequency (cycles/deg)
			d.freqT = (.6 * 2 .^ (0: 6))'; % cleaned temporal frequency (Hz)
			a = 10 .^ d.sensLog; d.amp = a; % contrast sensitivity (Hz/%)
		
			% Phase data and response
			folder = [m.folder, filesep, 'Yeh (95) phase'];
			dC = m.readFolder(folder); % response phase versus temporal frequency
			p = dC.phase; d.phase = p; % response phase (deg)
			d.resp = a .* exp(1i * pi * p / 180); % response (Hz), complex: fs x 1

			% Describe
			d.Properties.VariableDescriptions{'amp'} = 'Contrast sensitivity (Hz/%)';
			d.Properties.VariableDescriptions{'phase'} = 'Response phase (deg)';

	end
	d.Properties.VariableDescriptions{'freqS'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'freqT'} = 'Temporal frequency (Hz)';

function [d, m] = respGang(~, m) % ganglion cell freq. resp.: Wool (18)

	% Create table for L-cone stimulation
	f = [m.folder, filesep, 'Wool (18)/']; % data folder
	load([f, 'R1.mat'], 'R1'); % log amplitude (Hz)
	load([f, 'Rphase.mat'], 'Rphase'); % phase (deg)
	fs = size(R1, 1); % number of frequencies
	c = repmat([.45, 0], [fs, 1]); % contrast
	n = {'freqLogA', 'freqLogP', 'actLog', 'phase', 'cont'};
	dR = table(R1(:, 1), Rphase(:, 1), R1(:, 2), Rphase(:, 2), c, ...
		'variableNames', n);

	% Add data for M-cone stimulation
	load([f, 'G1.mat'], 'G1'); % log amplitude (Hz)
	load([f, 'Gphase.mat'], 'Gphase'); % phase (deg)
	fs = size(G1, 1); % number of frequencies
	c = repmat([0, .45], [fs, 1]); % contrast
	dG = table(G1(:, 1), Gphase(:, 1), G1(:, 2), Gphase(:, 2), c, ...
		'variableNames', n);

	% Process and store the combined data
	d = [dR; dG]; % combine tables
	d.freq = 10 .^ mean([d.freqLogA, d.freqLogP], 2); % spat. freq. (c. / deg)
	d.act = 10 .^ d.actLog; % impulse rate (Hz)
	d = sortrows(d, {'cont', 'freq'}); % sort on the independent variables
	d.Properties.VariableDescriptions = {'', '', '', ...
		'Phase (deg)', 'Contrast', 'Spatial frequency (cycles / deg)', ...
		'Impulse rate (Hz)'};

function [d, m] = sepBeta(d, m) % separation between 4CBeta cells: pragmatic

	%	Calculate density of 4CBeta cells from density of ganglion cells
	import prim.readTab % find function
	e = [0, 1, 3, 10, 30]; % eccentricities commonly used in runPrim (deg)
	if isfield(m.sepBeta, 'mult') % beta cell density is multiple of g.c. density

		%	Calculate densities
		dGangOff = readTab('densGangOff', e, m); % off-g.c. dens. (deg^-2): 1 x es
		dGangOn = readTab('densGangOn', e, m); % on-g.c. dens. (deg^-2): 1 x es
		dGang = dGangOff + dGangOn; % ganglion cell density (deg^-2): 1 x es
		dBeta = m.sepBeta.mult * dGang; % pragmatic beta dens. (deg^-2): 1 x es
		dEx = (1 - m.p.ratBetaIn) * dBeta; % excit. cells dens. (deg^-2): 1 x es
		dIn = m.p.ratBetaIn * dBeta; % dens. of inhibitory cells (deg^-2): 1 x es
	
		%	Calculate separation of 4CBeta cells and display
		sepEx = 1 ./ sqrt(dEx); % excit. sep., assuming square array (deg): 1 x es
		sepIn = 1 ./ sqrt(dIn); % inhibitory separation (deg): 1 x es

	else % fixed number of excitatory cells

		%	Calculate excitatory cell separation
		w = readTab('width', e); % visual field width (deg): 1 x es
		nE = m.sepBeta.num; % number of excitatory cells
		n = sqrt(nE); % number of cells on a visual field side
		sepEx = w / n; % excitatory cell separation (deg): 1 x es

		%	Calculate inhibitory cell separation
		r = m.p.ratBetaIn; % ratio of inhibitory to all beta neurons
		nI = (r / (1 - r)) * nE; % number of inhibitory cells
		n = sqrt(nI); % number of cells on a visual field side
		sepIn = w / n; % inhibitory cell separation (deg): 1 x es
		
	end
	fprintf('%6.4f, ', sepEx); fprintf('\n'); % print separations
	fprintf('%6.4f, ', sepIn); fprintf('\n'); % print separations

function [d, m] = showIm(~, m) % show natural images

	% Set image properties
	import prim.readIm prim.readImFolder;
	if isfield(m.showIm, 'files') % serial numbers of image files
		file = m.showIm.files; % user-specified
	else
		file = []; % default is all
	end
	[folder, name] = readImFolder(file); % list names of image files
	fs = length(name); % number of image files
	%	if isempty(file), file = 1: fs; end % serial numbers of all image files
	k = 0; % index of eligible file
	d = table; % empty table

	%	Loop over files
	for i = 1: fs % loop over files

		%	Read file and skip if ineligible
		nameC = name(i); % name of current file
		%	figure('windowStyle', 'docked'); imshow(folder + filesep + nameC);
		[im, sOrig] = readIm(folder, nameC); % read image file: ys x xs x cs
		if isempty(im) % file is unreadable, too small, or greyscale
			continue
		else
			k = k + 1; % row number
			s = size(im); % image size
			x = [1, s(2)]; y = [1, s(1)]; % limits for image display
			d(k, :) = table(x, y, sOrig, nameC, {im}, ...
				'variableNames', {'x', 'y', 'size', 'name', 'image'});
		end

	end

function [d, m] = wid(d, m) % find visual field width for fixed number of cones

	n = m.wid.cones; % required number of cones
	dens = d.densDeg; % cone density (cones/deg^2): es x 1
	d.width = sqrt(n ./ dens); % visual field width (deg): es x 1

function [d, m] = xyz2lms(d, m) % convert from XYZ colour space to LMS space

	switch 'stockman' % source
		case 'stockman' %	Stockman (23): 2° XYZ, Stockman (00) cone fundamentals
			l2x = [	1.94735469, -1.41445123,	0.36476327
							0.68990272,	0.34832189,		0
							0,          0,						1.93485343]; % lms2xyz
			x2l = inv(l2x); % invert matrix: xyz2lms
		case 'smith' % Smith (75) as interpreted by
			% daltonlens.org/understanding-cvd-simulation/#smith_spectral_1975
			x2l = [	.15514,		.54312, -.03286
							-.15514,	.45684,	.03286
							0,				0,			.01608];
	end
	m.p.xyz2lms = x2l; % insert this value into setLit.m

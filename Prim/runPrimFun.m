%	Calculate model responses

function m = runPrimFun % define function handles for functions in this file

	m.bandwidth = @bandwidth;
	m.calConv = @calConv;
	m.calStim = @calStim;
	m.calStruct = @calStruct;
	m.cor = @cor;
	m.dir2orient = @dir2orient;
	m.doCrop = @doCrop;
	m.getVal = @getVal;
	m.ifftReal = @ifftReal;
	m.listLoc = @listLoc;
	m.oi = @oi;
	m.osi = @osi;
	m.readConst = @readConst;
	m.readFile = @readFile;
	m.rect = @rect;
	m.rgb2contOld = @rgb2contOld;
	m.samp = @samp;
	m.saturate = @saturate;
	m.setStruct = @setStruct;
	m.setSyn = @setSyn;
	m.setVal = @setVal;
	m.setWeight = @setWeight;
	m.showOsi = @showOsi;
	m.solveF = @solveF;
	m.solveT = @solveT;

function loc = array(wid, sep, off) % calculate triangular array

%	Input:
%		wid = visual field width (deg)
%		sep = distance between nearest-neighbour nodes (deg)
%	Optional input:
%		off = offset, [x, y] (deg)
% Output:
%		loc = node locations, [x, y] (deg): ls x 2, where ls is the number of nodes
% Method:
%		The code calculates a parallelogram of nodes and then crops it to a square.
%		The parallelogram is wide enough to fill the visual field when cropped.

	%	Generate a parallelogram of nodes
	if ~ exist('off', 'var'), off = [0, 0]; end % default offset
	w = .5 * wid; % half-width (deg)
	u = exp(1i * pi / 3); % unit vector along oblique axis
	d = (1 + real(u)) * w + off(1); % distance to right of centre (deg)
	i = ceil(d / sep); % number of nodes to right of centre
	x = linspace(- i, i, 2 * i + 1) * sep; % hor. dist. of nodes from centre (deg)
	d = w + off(2); % distance above centre (deg)
	i = ceil(d / (imag(u) * sep)); % number of nodes above centre
	y = linspace(- i, i, 2 * i + 1) * sep; % ver. dist. of nodes from centre (deg)
	loc = x' + y * u; % generate array of nodes
	
	% Turn the nodes into locations and trim to a square
	loc = [real(loc(:)), imag(loc(:))]; % locations (deg): ls x 2, ls too big
	loc = loc + off; % include offset
	i = abs(loc); i = i(:, 1) <= w & i(:, 2) <= w; % nodes within square
	loc = loc(i, :); % locations of nodes within square (deg): ls x 2

function w = bandwidth(x, z, zCrit, cyclic) % calculate tuning bandwidth

%	Input:
%		x = stimulus: 1 x xs
%		z = response: xs x ls
%		zCrit = response level at which to calculate bandwidth: 1 x ls
%		cyclic = 1 for cyclic stimulus, 0 otherwise
% Assumption:
%		xs are equally spaced

	%	Loop over cell locations
	if cyclic % stimulus is cyclic: convert from closed to open interval
		z(end, :) = []; %	remove repeated point at end: xs -> xs - 1
		iCen = floor(.5 * size(z, 1)); % index of central point
	end
	ls = size(z, 2); % number of cell locations
	w = nan(2, ls); % allocate storage: left- and right-bandwidths: 2 x ls
	for loc = 1: ls % loop over cell locations

		%	Loop over half-bandwidths
		zC = z(:, loc); % current response: xs x 1
		[~, iMax] = max(zC); % current response maxiumum
		if cyclic % stimulus is cyclic: centre maximum within interval
			zC = circshift(zC, iCen - iMax); % centre maximum within interval
			iMax = iCen; % new index for maximum
		end
		for j = 1: 2 % left then right half-bandwidth
		
			%	Remove all but descending arm
			if j == 1 % left decending arm
				zD = zC(1: iMax); zD = flip(zD); % flip left for right
			else % right descending arm
				zD = zC(iMax: end);
			end

			% Interpolate across criterion response level
			zCritC = zCrit(loc); % current response criterion
			i = zD > zCritC; % indices of values above bandwidth criterion
			i = find(~ i, 1); % index of first value below criterion
			i = i - 1: i; % indices containing criterion: use 1: i?
			if ~ isempty(i) % avoid cases in which response doesn't fall below crit.
				i = interp1(zD(i), i, zCritC); % criterion index
				w(j, loc) = (i - 1) * (x(2) - x(1)); % bandwidth
			end

		end
		
	end

function m = betaEx(m) % calculate 4CBeta excitatory array

	% Calculate separation between neighbouring cells
	import prim.readTab % find function
	switch 'gang' % choose method
		case 'cort' % set separation from cortical density
			sep = m.p.densBeta; % 4CBeta neuronal density (mm^-2): 1 x 1
			sep = (1 / m.p.ratCort) ^ 2 * sep; % density (deg^-2)
			sep = m.p.kDensBeta * sep; % reduce density to shorten computation time
			sep = (1 - m.p.ratBetaIn) * sep; % density of excitatory cells (deg^-2)
			sep = 1 / sqrt(sep); % separation, assuming square array (deg)
		case 'gang' % set separation from ganglion cell density
			% Justification: Merigan (90) showed macaque psychophysical resolution was
			% predictable from ganglion cell density of either sign
			sep = readTab('sepEx', m.p.ecc); % separation of excitatory cells (deg)
	end

	% Calculate cell locations
	x = 0: sep: .5 * m.p.wid; % non-negative x values (deg)
	xN = - fliplr(x(2: end)); % negative x values (deg)
	x = [xN, x]; % x values (deg)
	[x, y] = ndgrid(x); % x and y values for all array locations (deg)
	m.p.betaEx.loc = [x(:), y(:)]; % cell locations (x, y) (deg): ls x 2
	m.p.betaEx.type = "ex"; % cortical cell type: 1 x 1

function m = betaIn(m) % calculate 4CBeta inhibitory array and convergence fun.

	% Calculate separation between neighbouring cells
	import prim.readTab % find function
	switch 'gang' % choose method
		case 'cort' % set separation from cortical density
			sep = m.p.densBeta; % 4CBeta neuronal density (mm^-2): 1 x 1
			sep = (1 / m.p.ratCort) ^ 2 * sep; % density (deg^-2)
			sep = m.p.kDensBeta * sep; % reduce density to shorten computation time
			sep = m.p.ratBetaIn * sep; % density of inhibitory cells (deg^-2)
			sep = 1 / sqrt(sep); % separation between neighbouring cells (deg)
		case 'gang' % set separation from ganglion cell density
			sep = readTab('sepIn', m.p.ecc); % separation of excitatory cells (deg)
	end

	% Calculate cell locations
	xs = floor(.5 * m.p.wid / sep); % number of positive x values
	widH = xs * sep; % half width of grid (deg)
	x = widH * linspace(-1, 1, 2 * xs + 1); % x values (deg)
	[x, y] = ndgrid(x); % x and y values for all array locations (deg)
	m.p.betaIn.loc = [x(:), y(:)]; % cell locations (x, y) (deg): ls x 2
	m.p.betaIn.type = "in"; % cortical cell type: 1 x 1

function g = calConv(locI, locO, rad) % calculate convergence matrix

%	Input:
%		locI = locations (x, y) of input neurons (deg): is x 2
%		locO = locations (x, y) of output neurons (deg): os x 2
%		rad = radius of Gaussian attenuation (deg): 1 x 1
%	Output:
%		g = Gaussian function of distance between locI and locO (deg): os x is

	locI = permute(locI, [3, 1, 2]); % prepare for subtraction: 1 x is x 2
	locO = permute(locO, [1, 3, 2]); % os x 1 x 2
	g = locO - locI; % displacement of output neuron from input: os x is x 2
	g = sum(g .^ 2, 3); % squared distance of output neuron from input: os x is
	if rad > 0 % standard case
		g = exp(- g / rad ^ 2); % Gaussian function of distance: os x is
	else % private line
		g = eye(size(g, 1)); % Gaussian is an impulse: os x is
	end

function dpot = calDeriv(t, pot, m) % calc. derivatives of model equations

	% Sort potential into cell types
	i = 0; % current row in potential vector
	for c = m.p.cell % loop over cell types
		type = c.type; % current type
		a = c.array; % location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(type) = pot(i + (1: ls)); % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end

	% Calculate derivatives for cone array
	rate = 1 / m.p.tau; % reciprocal of time constant (s^-1)
	drive = driveT(t, m)'; % cone generator potential (mV): ls x 1
	back = m.p.kSur * m.p.wHorCone * p.hor; % feedback potential (mV): ls x 1
	if m.p.coneStages == 1 % one cone stage
		dp.cone = rate * (- m.p.kSens * drive - m.p.back * back - p.cone);
			% cone potential (mV/s): ls x 1
	else % two cone stages
		dp.coneInt = rate * (- m.p.kSens * drive - m.p.back * back - p.coneInt);
			% intermediate cone potential (mV/s): ls x 1
		dp.cone = rate * (p.coneInt - p.cone); % cone potential (mV/s): ls x 1
	end
	dp.hor = rate * (m.p.wConeHor * p.cone - p.hor); % hor. cell (mV/s): ls x 1
	dp.back = 0 * dp.hor; % feedback (mV/s): ls x 1
	dp.biOff = rate * (p.cone - p.biOff); % off-bipolar (mV/s): ls x 1
	dp.biOn = rate * (- p.cone - p.biOn); % on-bipolar (mV/s): ls x 1
	
	% Calculate derivatives for ganglion cell arrays
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i; gen = "gen" + i;
			% construct names
		dp.(gang) = rate * (m.p.(w) * p.(bi) + m.p.potRest - p.(gang));
			% ganglion cell potential (mV): ls x 1
		pC = max(0, p.(gang)); % presynaptic pot. at geniculate (mV)
		dp.(gen) = rate * (m.p.kGangGen * pC - p.(gen)); % gen. pot. (mV): ls x 1
	end

	% Calculate derivatives for layer 4CBeta
	pGen = [p.genOff; p.genOn]; % combine off- and on-inputs: ls x 1
	pGen = max(0, pGen); % presynaptic potential, rectified by the geniculate
	pEx = m.p.kGenBeta * m.p.w_gen_betaEx * pGen; % weighted and summed exc. input
	pIn = m.p.kGenBeta * m.p.w_gen_betaIn * pGen;
	dp.betaIn = rate * (pIn - p.betaIn); % deriv. of inhib. pot.: is x 1
	pIn = max(0, p.betaIn); % inhibitory potential, rectified by inh. cells
	pIn = m.p.kInEx * m.p.w_betaIn_betaEx * pIn; % weighted, summed inhib. input
	dp.betaEx = rate * (pEx - pIn - p.betaEx); % deriv. of ex. pot. (mV/s): ls x 1
	dp.sumGen = 0 * dp.betaEx; dp.sumIn = 0 * dp.betaEx; % place-holders

	%	Vectorise derivatives
	dpot = zeros(size(pot)); % initialise derivatives
	i = 0; % current row in potential vector
	for c = m.p.cell % loop over cell types
		type = c.type; % current type
		dpC = dp.(type); ls = size(dpC, 1); % derivatives, number of derivatives
		dpot(i + (1: ls)) = dpC; % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end

function s = calStim(t, loc, m) % calculate stimulus movie

% Inputs:
%		t: times (s): vector with length ts
%		loc = visual field locations, (x, y) (deg): ls x 2
%		m = metadata
%	Ouput:
%		s(t, loc, c) (contrast-units): ts x ls x cs

	% Initialise
	dir = (pi / 180) * m.p.dir; % motion dir. (radians ACW from rightward)
	fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
	fT = 2 * pi * m.p.freqT; % temporal frequency (radians/s)
	pS = (pi / 180) * m.p.phaseS; % spatial phase (radians)

	% Transform location to distance in direction of motion
	u = loc * [cos(dir); sin(dir)]; % dist. relative to centre (deg): ls x 1
	u = u'; % make it a row vector (deg): 1 x ls
	u = u - m.p.locS; % distance relative to stimulus location (deg): 1 x ls

	% Calculate stimulus
	t = t(:); % make time a column vector: ts x 1
	switch m.p.stimT % drift or pulse
		case 'drift' % drifting grating
			s = cos(fS * u - pS - fT * t); % stimulus: ts x ls
		case 'pulse' % pulsed grating
			s = (t <= m.p.dur) .* cos(fS * u - pS); % stimulus: ts x ls
	end

	% Multiply by contrast
	c = m.p.cont; % contrast: 1 x cs
	c = shiftdim(c, -1); % prepare to multiply: 1 x 1 x cs
	s = c .* s; % multiply: ts x ls x cs

function m = calStruct(m) % calculate the model structure

	% Initialise randomisation and warnings
	rng(m.p.seed); % make the results reproducible
	i = 'runPrim:LargeImaginary'; % warning about large imag. components in ifft
	if m.p.warning, warning('on', i); % define warning, turn it on
	else, warning('off', i); end % define warning, turn it off
	
	%	Create a list of stages
	a = m.p.array; % cell array of neuronal array names
	c = cell(size(a)); % cell array to hold structure for each neuronal array
	for i = 1: length(a) % loop over neuronal arrays
		aC = a{i}; % name of current array
		c{i} = struct('type', m.p.(aC).stage, 'array', aC); % stuct. for this array
	end
	m.p.cell = horzcat(c{:}); % concatenate across structures

	% Calculate the subcortical convergence function radii
	import prim.readTab; % function location
	e = m.p.ecc; % current eccentricity (deg)
	m.p.radOpt = polyval(m.p.radOptCoef, e); % point spread function (deg)
	m.p.radHor = exp(polyval(m.p.radHorCoef, e)); % hor. cell rec. field (deg)
	%	m.p.radGang = polyval(m.p.radGangCoef, e); % ganglion cell dendritic (deg)
	m.p.radGang = readTab('radGang', e); % ganglion cell dendritic (deg)
	m.p.radCen = sqrt(m.p.radOpt ^ 2 + m.p.radGang ^ 2); % centre radius (deg)
	m.p.radSur = sqrt(m.p.radOpt ^ 2 + 2 * m.p.radHor ^ 2); % surround rad. (deg)

	% Calculate the cortical convergence function radii
	fac = readTab('ratCort', e); % inverse cort. mag'n factor (deg/mm)
	m.p.ratCort = fac; % store
	m.p.radBeta = fac * m.p.radBetaMm; % 4CBeta radius (deg)
	m.p.radInEx = fac * m.p.radInExMm; % inhibitory-excitatory radius (deg)
	m.p.radBetaCen = sqrt(m.p.radCen ^ 2 + m.p.radBeta ^ 2); % beta centre (deg)
	m.p.radBetaSur = sqrt(m.p.radSur ^ 2 + m.p.radBeta ^ 2); % beta surround (deg)

	% Calculate the arrays and convergence functions
	m = cone(m); % cone array
	m = gang(m); % ganglion cell array
	m = betaEx(m); % 4CBeta excitatory array
	m = betaIn(m); % 4CBeta inhibitory array
	m = weight(m); % convergence functions

	%	Calculate the subcortical resting potential
	m.p.potRest = m.p.actRest / m.p.kRect; % ganglion cell resting potential (mV)

function m = cone(m) % calculate the cone array

	% Calculate the cone locations
	import prim.readTab; % function location
	dens = readTab('densCone', m.p.ecc); % density at specified ecc. (deg^-2)
	s = 1 / sqrt(sin(pi / 3) * dens); % cone separation (deg)
	loc = array(m.p.wid, s); % calculate the cone locations (x, y) (deg): cs x 2
	m.p.cone.loc = loc; % cone locations (x, y) (deg): cs x 2
		
	% Assign L and M cone types to all locations, ignoring S cones
	cs = size(loc, 1); % number of cones
	type = ones(cs, 1); % set all cones to L type: cs x 1
	ts = cs * m.p.ratCone; % number of L, M, S cones
	ls = ts(1); ms = ts(2); ss = ts(3); % number of cones of each type
	ms = round((ms / (ls + ms)) * cs); % number of M cones
	i = randperm(cs); % randomise cone order: 1 x cs
	type(i(1: ms)) = 2; % set fraction of cones to M type: cs x 1
	
	% Make an offset triangular array for the S cones and reassign these loc. to S
	dens = (ss / cs) * dens; % S cone density (cells/deg^2)
	s = 1 / sqrt(sin(pi / 3) * dens); % S cone separation (deg)
	off = .5 * s * tan(pi / 6) * [1, 1]; % S cone array offset, [x, y] (deg)
	locS = array(m.p.wid, s, off); % S cone array
	i = knnsearch(loc, locS); % find the nearest locations in the L, M array
	type(i) = 3; % set those locations to S type: cs x 1
	m.p.cone.type = type; % cone types: cs x 1

function c = cor(m, pIn, pOut) % cross-correlate synaptic inputs and outputs

%	Inputs: response time course in Fourier domain (mV):
%		pIn = input: ts x js
%		pOut = output: ts x ks

	pIn = permute(pIn, [1, 3, 2]); % prepare for cross-correlation: ts x 1 x js
	c = conj(pIn) .* pOut; % transform of cross-correlation (mV^2): ts x ks x js
	c = ifft(c); % transform to time domain (mV^2): ts x ks x js
	lagWin = m.p.t >= 0 & m.p.t <= .1; % cross-correlation lags from 0 to .1 s
	c = c(lagWin, :, :); % keep lags from 0 to .1 s (mV^2): ts x ks x js, new ts
	[~, lag] = max(abs(c), [], 'linear'); % lag at extremes: 1 x ks x js
	c = c(lag); % extreme cross-correlation (mV^2): 1 x ks x js
	c = shiftdim(c, 1); % remove single-element dimension: ks x js

function orient = dir2orient(dir) % convert stimulus direction to cyclic orient.

%	Input: dir = stimulus direction (deg), [-180, 180): vector
%	Output: orient = complex vector: length = 1, angle = 2 * dir

	o = pi * dir / 180; % convert to radians
	o = 2 * o; % convert to cyclic
	orient = exp(1i * o); % make it complex

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

function d = driveF(m) % calculate Fourier transform of cone drive

%	Output:
%		d = Fourier transform of drive for each cone (mV): ls x fs
%	Method:
%		Calculate stimulus as a function of visual field location
%		Add temporal frequency component
%		Cross-correlate stimulus and optical point spread function
%		Set contrast

	%	Initialise
	import prim.transIm % find function
	type = m.p.cone.type; % cone type: ls x 1

	% Calculate stimulus
	switch m.p.stimS % stimulus type
		case 'grating' % grating

			% Initialise
			fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
			pS = (pi / 180) * m.p.phaseS; % spatial phase (radians)
			loc = m.p.cone.loc; % cone location (deg): ls x 2
			ls = size(loc, 1); % number of cones
			f = 2 * pi * m.p.f; % fft frequencies (radians/s): 1 x fs
			fs = length(f); % number of frequencies

			%	Calculate distance
			dir = (pi / 180) * m.p.dir; % motion dir. (radians ACW from rightward)
			u = loc * [cos(dir); sin(dir)]; % dist. relative to centre (deg): ls x 1
			u = u - m.p.locS; % distance relative to stimulus location (deg): ls x 1

			%	Calculate Fourier transform of stimulus
			switch m.p.stimT % stimulus temporal component
				case 'drift' % drifting grating
					s = exp(- 1i * (fS * u - pS)); % phase shift at cone location: ls x 1
					d = zeros(ls, fs); % cone signal, all frequencies: ls x fs
					d(:, 2) = s; % fundamental component (mV): ls x fs
					d(:, end) = conj(s); % negative fundamental (mV): ls x fs
					d = .5 * fs * d; % fundamental fft component (mV): ls x fs
				case 'pulse' % pulsed grating
					s = cos(fS * u - pS); % stimulus: ls x 1
					d = s .* pulse(m); % multiply by pulse transform: ls x fs
			end

			%	Include optical spread and cone contrast
			d = exp(-.25 * (m.p.radOpt * fS) ^ 2) .* d; % optics: ls x fs
			%{
			for i = 1: ls % loop over cones
				j = m.p.cone.type(i); % cone type
				d(i, :) = m.p.cont(j) * d(i, :); % multiply by appropriate cone contrast
			end
			%}
			d = m.p.cont(type)' .* d; % multiply by appropriate cone contrast

		case 'image' % image
			d = transIm(m); % transduce image to cone drive (mV): ls x ts
	end

	%	Adapt L-cones
	if m.p.adapt % adapt L-cones
		i = type == 1; % indices of L-cones: ls x 1
		d(i, :) = m.p.kAdapt * d(i, :); % adapt them
	end

function d = driveT(t, m) % calculate cone drive at specified time(s)

	% Calculate stimulus and attenuate with the optical point spread function
	switch m.p.stimS % grating or image
		case 'grating'
			loc = m.p.cone.loc; % cone location (deg): ls x 2
			d = calStim(t, loc, m); % stimulus (contrast units): ts x ls x cs
			fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
			d = exp(-.25 * (m.p.radOpt * fS) ^ 2) .* d; % attenuated: ts x ls x cs			
		case 'image'

			% Image: read, convert to cone contrast, blur, and sample
			d = interpIm(m); % cone drive: ls x 1 x cs
			
			%	Add temporal component
			d = reshape(d, 1, [], 3); % prepare for multiplication: 1 x ls x cs
			d = (t <= m.p.dur) .* d; % multiply by temporal waveform: ts x ls x cs
			
	end

	% Choose cone contrast appropriate to each one
	[~, ls, cs] = size(d); % number of locations, cone types
	type = m.p.cone.type; % cone type: ls x 1
	i = sub2ind([ls, cs], (1: ls)', type); % convert to indices: ls x 1
	d = d(:, i); % cone-specific stimulus: ts x ls
	if m.p.adapt % adapt L-cones
		i = type == 1; % indices of L-cones
		d(:, i) = m.p.kAdapt * d(:, i); % adapt them: ts x ls
	end

function m = gang(m) % calculate the ganglion cell arrays

	import prim.readTab; % function location
	for s = ["Off", "On"] % off- then on-centre
		name = "densGang" + s; % density to read from table
		densC = readTab(name, m.p.ecc); % density at specified ecc. (deg^-2)
		sep = sqrt(1 / (sin(pi / 3) * densC)); % cell separation (deg)
		loc = array(m.p.wid, sep); % calculate cell locations (x, y) (deg): ls x 2
		switch m.p.align % adjust locations
			case 1 % align ganglion cells with bipolar cells
				locCone = m.p.cone.loc; % cone locations (deg): cones x 2
				j = knnsearch(locCone, loc); % indices of nearest cones: cs x 1
				loc = locCone(j, :); % match location to that of nearest cone: cs x 2
			otherwise % randomise ganglion cell locations
				dev = readTab('devGang', m.p.ecc, m); % stand. dev. about tri. array
				loc = loc + normrnd(0, dev * sep, size(loc)); % perturb (deg)
		end
		n = "gang" + s; % name of current array
		m.p.(n).loc = loc; % ganglion cell location (x, y) (deg): ls x 2
		m.p.(n).type = lower(s); % ganglion cell type: 1 x 1
	end

function [name, val, vals] = getVal(m, task) % get values of the stimulus par.

%	Outputs:
%		name = names of specified stimulus parameters: cell
%		val = values for the specified stimulus parameters:
%			cell array with one member for each stimulus parameter
%		vals = number of values for each specified stimulus parameter: double

	import prim.getIm; % find function
	switch m.p.stimS % grating or image
		case 'grating'

			% Find the stimulus parameters to be varied
			name = fieldnames(m.(task)); % names in m.task
			nameAll = fieldnames(m.p); % names of stimulus parameters
			i = ismember(name, nameAll); % m.task fields that are stimulus parameters
			name = name(i)'; % names of specified parameters: 1 x ns
			names = length(name); % number of parameters
			
			% Find the values of the response parameters, which may be multi-column
			val = cell(1, names); % values of parameters
			vals = zeros(1, names); % number of values
			for i = 1: names % loop over response parameters
				nameC = name{i}; % name of current parameter
				valC = m.(task).(nameC); % values of parameter
				val{i} = valC; % store values
				vals(i) = size(valC, 1); % update number of values
			end

		case 'image'
			[name, val, vals] = getIm(m, task); % get images to be used as stimuli
	end

function r = ifftReal(r, dim) % calculate real part of inverse transform

	if nargin == 1 % the frequency dimension is unspecified
		dim = 1; % assume that the frequency dimension is 1
	end
	r = ifft(r, [], dim); % take inverse transform
	s = warning('query', 'runPrim:LargeImaginary'); % check state of warning
	if s.state == "on" % check whether warning is required
		mReal = max(abs(real(r(:)))); % maximum of real components
		mImag = max(abs(imag(r(:)))); % maximum of imaginary components
		ratio = mImag / mReal; % ratio of imaginary to real maxima
		if ratio > 1e-5 % issue warning
			i = 'runPrim:LargeImaginary'; % warning ID
			w = 'Imaginary components of inverse transform are relatively large: %g';
			warning(i, w, ratio); % issue warning
			warning('off', i); % turn off further warnings with this ID
		end
	end
	r = real(r); % make it real

function loc = listLoc(m) % create a list of visual field locations

	[x, y] = ndgrid(m.p.x); % grid locations (deg): xs x ys
	loc = [x(:), y(:)]; % locations (x, y) (deg): ls x 2

function ind = oi(fL, fM) % calculate opponency index

	dif = fL - fM; % difference of fields: ls x ps
	indMin = min(dif, [], 2); indMax = max(dif, [], 2); % bounds of dif.: ls x 1
	num = indMax + indMin; den = indMax - indMin; % sum and dif. of bounds
	ind = 1 - abs(num ./ den); % index: low for single opponency, high for double

function ind = osi(dir, r) % calculate orientation selectivity index

%	Input:
%		dir = stimulus direction (deg), [-180, 180): 1 x ds
%		r = direction response (mV or Hz): ls x ds
%	Output:
%		ind = orientation selectivity index = 1 - circular variance, [0, 1]: ls x 1
%		oPref = preferred orientation (deg), [-90, 90): ls x 1 *** add? ***

	o = dir2orient(dir); % convert direction to orientation vector: 1 x ds
	ind = sum(r .* o, 2) ./ sum(r, 2); % o.s.i., complex: ls x 1
	ind = abs(ind); % o.s.i.: ls x 1

function P = pulse(m) % calculate pulse transform

	dur = m.p.dur; % pulse duration (s)
	f = 2 * pi * m.p.f; % fft frequencies (radians/s): 1 x fs
	fs = length(f); % number of frequencies
	a = (dur / 2) * f; % argument for pulse transform (rad): 1 x fs
	P = 2 * sin(a) ./ f; % pulse transform (s): 1 x fs
	P = exp(- 1i * a) .* P; % shift so that pulse starts at 0 (s): 1 x fs
	P(1) = dur; % replace nan (s)
	P = (1 / m.p.time) * P; % normalise for stim. time (no unit): 1 x fs
	P = fs * P; % change to fft units: 1 x fs

function [d, m] = readFile(m) % set data folder, read or create data file

	% Set names of data folder and data file
	f = string(userpath) + filesep + 'Data' + filesep + m.p.project; % folder
	if exist(f, 'dir') % folder exists
		m.folder = f; % set it
	else % no such folder
		m.folder = string(cd); % default is current folder
	end
	file = m.folder + filesep + m.p.file; % file name
	
	% Read or create file
	if ~ isempty(file) && exist(file, 'file') % data file exists
		load(file, 'd'); % load data file
		n = d.Properties.VariableNames; % variable name, cell array
		if ismember('ecc', n), m.p.ecc = d.ecc(1); end % eccentricity (deg)
		if ismember('width', n), m.p.wid = d.width(1); end % v.f. patch width (deg)
	else
		d = table; % default is empty data file
	end

function r = rect(r, dim) % rectify a frequency-domain response

	if nargin == 1 % the frequency dimension is unspecified
		dim = 1; % assume that the frequency dimension is 1
	end
	rC = ifftReal(r, dim); % transform to time domain and make it real
	if min(rC, [], 'all') < 0 % time domain response has negative elements
		rC = max(0, rC); % rectify
		r = fft(rC, [], dim); % return to frequency domain
	end

function im = rgb2contOld(im, m) % convert sRGB to cone contrast *** remove ***

	% Convert sRGB to LMS
	im = rgb2xyz(im); % convert from sRGB to XYZ: double
	[xs, ys, cs] = size(im); % size of image
	im = reshape(im, [], 3); % prepare for multiplication: ls x 3
	im = im * m.p.xyz2lms'; % convert XYZ to LMS % (ls x cs) * (cs x cs) = ls x cs
	im = reshape(im, [xs, ys, cs]); % restore shape: xs x ys x 3
	
	%	Convert LMS to cone contrast
	sig = m.p.radFix / sqrt(2); % convert rad. to st. dev. (deg)
	sig = (sig / m.p.wid) * xs; % background radius (pixels)
	b = imgaussfilt(im, sig); % local background: xs x ys x cs
	im = (im - b) ./ b; % map (contrast units): xs x ys x cs, double

function m = samp(m) % calculate sample locations, times, and temporal freq.

	%	Calculate the sample locations
	m.p.x = .5 * m.p.wid * linspace(-1, 1, m.p.xs); % x values (deg): 1 x xs
	
	% Calculate the number of samples and fundamental frequency
	n = m.p.ts; % number of sample times and temporal frequencies
	nH = floor(.5 * n); % half the number, converted to integer
	n = 2 * nH; % make sure that the number of samples is even
	fund = 1 / m.p.time; % fundamental frequency (Hz)
	%	switch m.p.stimT, case 'drift', fund = m.p.freqT; end

	% Calculate the sample frequencies
	ind = linspace(- nH, nH, n + 1); % sample indices
	ind(end) = []; % make the sequence open-ended
	f = fund * ind; % sample frequencies (Hz): n x 1
	m.p.f = fftshift(f); % begin with zero frequency, for compatibility with fft

	% Calculate the sample times
	ind = linspace(0, 1, n + 1); % sample indices
	ind(end) = []; % make the sequence open-ended
	t = (1 / fund) * ind; % sample times (s): n x 1
	m.p.t = t; % begin with time 0, for intuitive use

function z = saturate(z, dir, dirTun, quant) % show OSI by saturation in plot

%	Inputs:
%		z = plotted value (preferred direction or maximum response): ls x 1
%		dir = stimulus direction (deg): 1 x ds
%		dirTun = direction tuning (mV or Hz): ls x ds
%		quant = 'dir' for preferred direction, or 'max' for maximum response
%	Ouput:
%		z = RGB representation of plotted value: ls x 3

	%	Calculate saturation
	s = dirTun; % direction tuning (mV or Hz): ls x ds
	s = osi(dir, s); % saturation = orientation selectivity index: ls x 1
	s = s / prctile(s, 95); % normalised saturation: ls x 1
	s = min(1, s); % limit to 1, [0, 1]: ls x 1

	switch quant % calculate RGB array
		case 'dir' % preferred direction
			%	calculate hue and assemble HSV array
			h = z(:); % preferred direction (deg), [-180, 180): ls x 1
			h = mod(h + 90, 180); % orientation (deg), [0, 180): ls x 1
				% addition of 90 deg aligns RGB and HSV representations
			h = h / 180; % hue, [0, 1): ls x 1
			z = ones(size(h)); % value = 1: ls x 1
			z = [h, s, z]; % HSV array: ls x 3
		case 'max' % maximum
			%	convert colour map to HSV and insert saturation, value
			z = z(:); % response maximum (mV or Hz): ls x 1
			c = colormap; % current colour map: cs x 3
			i = round((z / max(z)) * size(c, 1)); % index into map: ls x 1
			z = c(i, :); % read colour from colour map: ls x 3
			z = rgb2hsv(z); % convert to HSV: ls x 3
			z(:, 2) = s; % insert saturations: ls x 3
			z(:, 3) = 1; % insert values: ls x 3
	end
	z = hsv2rgb(z); % RGB array: ls x 3

function m = setSyn(m) % create list of all plastic synapses

	stage = [m.p.synHebb; m.p.synIn]; % one row for each site
	ss = size(stage, 1); % number of sites
	s = struct; % empty structure
	for i = 1: ss % loop over sites
		s(i).stage = stage(i, :); % stages in current site
		s(i).mech = "Hebb"; % mechanism of plasticity
		n = stage(i, 1) + "_" + stage(i, 2); % site name
		s(i).site = n; % store it
		s(i).att = "a_" + n; % name of distance-based attenuation
		s(i).weight = "w_" + n; % name of synaptic weight
		if i < ss % Hebbian synapse
			sC = size(m.p.(s(i).att)); % number of output and input neurons
			s(i).mod = ones(sC); % modulation factor: os x is
		else % last row is for the spike-timing-independent synapse
			s(i).mech = "In"; % in for independence and inhibition
			s(i).mod = 1; % % inhibitory-excitatory gain
		end
	end
	m.p.syn = s; % store structure: 1 x ss

function m = setVal(m, name, val, iVal) % set the stimulus parameter values

%	Inputs:
%		m = metadata
%		name = names of specified stimulus parameters
%		val = values for the specified stimulus parameters: cell array with
%			one member for each stimulus parameter
%		iVal = linear index of current set of values

	%	Determine the subscripts of the stimulus parameters
	import prim.setIm % find function
	names = length(name); % number of parameters
	switch names % calculate size of val
		case 0, return % nothing to do
		case 1, s = [length(val{1}), 1]; % vector
		otherwise, s = arrayfun(@(x) size(x{1}, 1), val);
			% number of values for each parameter
	end
	sub = cell(1, names); % value subscripts
	[sub{:}] = ind2sub(s, iVal); % subscripts for each parameter

	%	Set the values
	switch m.p.stimS % grating or image
		case 'grating'
			for i = 1: names % loop over parameters
				nameC = name{i}; % parameter's name
				valC = val{i}; % parameter's values
				subC = sub{i}; % subscript of the current value
				m.p.(nameC) = valC(subC, :); % set value for this stimulus parameter
			end
		case 'image'
			m = setIm(m, val, sub); % set up the image stimulus
	end

function m = setWeight(d, m) % use data in table d to set synaptic weights in m

	% Set the plastic synaptic weights
	if isempty(d) % data file empty: useful for subcortex without cortex
		m.p.w_gen_betaEx = m.p.a_gen_betaEx; % modulation = 1: ks x js
		m.p.w_gen_betaIn = m.p.a_gen_betaIn; % modulation = 1: ks x js
		m.p.kInEx = 1; % modulation = 1: 1 x 1
	else % data file contains modulations
		syn = d.syn; % synaptic sites
		syns = length(syn) - 1; % number of Hebbian sites
		for i = 1: syns % loop over synaptic sites
			att = syn(i).att; % name of distance-based attenuation
			w = syn(i).mod .* m.p.(att); % synaptic weights: ks x js
			w = w ./ sum(w, 2); % normalise to unity sum: ks x js
			m.p.(syn(i).weight) = w; % store: ks x js
		end
		m.p.kInEx = syn(end).mod; % inhibitory-excitatory gain: 1 x 1
	end

function p = solveF(m) % solve the model equations in the frequency domain

	% Initialise
	f = 2 * pi * m.p.f; % fft frequencies (radians/s): 1 x fs
	T = 1 ./ (1 + 1i * m.p.tau * f); % temporal attenuation (radians): 1 x fs
	D = - m.p.kSens * driveF(m); % drive (mV): ls x fs
	w = m.p.back * m.p.kSur * m.p.wHorCone * m.p.wConeHor; % spat. filter: ls x ls
	ls = size(w, 1); % number of cones
	n = m.p.coneStages; % number of cone stages
	fs = length(f); % number of frequencies

	%	Initialise for for-loop
	if any(D(:, 3: end - 1), 'all') % speed up calculation for drifting grating
		range = 1: fs; % all frequencies
	else
		range = [1, 2, fs]; % mean and fundamental components
	end
	pC = zeros(ls, fs); % allocate storage: ls x fs

	% Calculate cone signal	
	p.coneInt = pC; % intermediate cone signal; no need to calculate
	for j = range % loop over frequencies
		TC = T(j); % current temporal filter (radians): 1 x 1
		DC = D(:, j); % current drive (mV): ls x 1
		pC(:, j) = (eye(ls) + TC .^ (n + 1) * w) \ (TC .^ n .* DC);
			% cone signal (mV): ls x 1
	end
	p.cone = pC; % cone signal (mV): ls x fs
	
	% Calculate remaining cone array signals
	p.hor = T .* (m.p.wConeHor * p.cone); % horizontal cell signal (mV): ls x fs
	p.back = m.p.kSur * m.p.wHorCone * p.hor; % feedback signal (mV): ls x fs
	p.biOff = T .* p.cone; % off-bipolar signal (mV): ls x fs
	p.biOn = - T .* p.cone; % on-bipolar signal (mV): ls x fs
	
	% Calculate ganglion cell array signals
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i; gen = "gen" + i;
			% construct weight name
		p.(gang) = T .* (m.p.(w) * p.(bi)); % ganglion cell fund. (mV): ls x fs
		p.(gang)(:, 1) = p.(gang)(:, 1) + fs * m.p.potRest; % rest. pot.: ls x fs
		p.(gen) = m.p.kGangGen * T .* rect(p.(gang), 2); % gen. pot. (mV): ls x fs
	end

	% Calculate layer 4CBeta array signals
	pGen = [p.genOff; p.genOn]; % combine off- and on-centre inputs (mV): ls x fs
	pGen = rect(pGen, 2); % rectify input (mV): ls x fs
	pEx = m.p.kGenBeta * m.p.w_gen_betaEx * pGen; % summed exc. input (mV): ls x fs
	p.sumGen = pEx; % store summed excitation (mV): ls x fs
	pIn = m.p.kGenBeta * m.p.w_gen_betaIn * pGen;
	p.betaIn = T .* pIn; % inhibitory cell response (mV): ls x fs
	pIn = m.p.kInEx * m.p.w_betaIn_betaEx * rect(p.betaIn, 2); % summed: ls x fs
	p.sumIn = pIn; % store summed inhibition (mV): ls x fs
	p.betaEx = T .* (pEx - pIn); % excitatory cell response (mV): ls x fs

	% Return stages, not arrays
	s = m.p.cell; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).type; % current stage
		s(i).resp = p.(sC); % add response
	end
	p = s; % 1 x zs

function p = solveT(m) % solve the model equations in the time domain

	% Set initial values for potentials
	for c = m.p.cell % loop over stages
		s = c.type; % current stage
		a = c.array; % current location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(s) = zeros(ls, 1); % default resting potential (mV): ls x 1
	end
	pRest = m.p.potRest; % ganglion cell resting potential (mV)
	p.gangOff(:) = pRest; p.gangOn(:) = pRest; % ganglion cells
	pRest = m.p.kGangGen * pRest; % geniculate resting potential (mV)
	p.genOff(:) = pRest; p.genOn(:) = pRest; % geniculate
	pRest = m.p.kGenBeta * pRest; % 4CBeta resting potential (mV)
	p.betaIn(:) = pRest; % 4CBeta inhibitory cells
	p.betaEx(:) = (1 - m.p.kInEx) * pRest; % 4CBeta excitatory cells
	pI = struct2cell(p); % make it a cell array
	pI = vertcat(pI{:}); % concatenate across stages: ls x 1, all locations

	% Numerically integrate the equations
	fun = @(t, p)calDeriv(t, p, m); % function to calculate derivatives
	t = m.p.t; % sampling times (s): 1 x ts
	[~, pot] = ode45(fun, t, pI); % potential (mV): ts x ls
	
	% Unpack potential
	i = 0; % current column in potential
	for s = m.p.cell % loop over stages
		sC = s.type; % current stage
		ls = length(p.(sC)); % number of locations
		pC = pot(:, i + (1: ls)); % store potentials in named array: ts x ls
		p.(sC) = pC'; % transpose to standard dimension order: ls x ts
		i = ls + i; % update column
	end

	% Add unintegrated signals
	p.back = m.p.kSur * m.p.wHorCone * p.hor; % f'back (mV):
		%	(ls x ls) x (ls x ts) = ls x ts
	pExc = [p.genOff; p.genOn]; % combine off- and on-inputs
	pExc = max(0, pExc); % presynaptic potential, rectified by the geniculate
	p.sumGen = m.p.kGenBeta * m.p.w_gen_betaEx * pExc; % summed excit. input
	pInh = max(0, p.betaIn); % inhibitory potential, rectified by inh. cells
	p.sumIn = m.p.kInEx * m.p.w_betaIn_betaEx * pInh; % summed inhibitory input

	% Rearrange potential to match old format
	s = m.p.cell; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).type; % current stage
		s(i).resp = p.(sC); % add response: ls x ts
	end
	p = s; % structure, including responses: 1 x zs

function m = weight(m) % calculate convergence function synaptic weights

	% Calculate cone to horizontal attenuation
	locI = m.p.cone.loc; % cone locations: cs x 2
	locO = locI; % horizontal cells are colocated with cones: cs x 2
	r = m.p.radHor; % horizontal cell field radius: 1 x 1
	a = calConv(locI, locO, r); % distance-based attenuation: cs x cs
	
	% Set the weight to zero for absent inputs, and store
	iCone = m.p.cone.type == 3; % cone is S-type: cs x 1
	w = a; w(:, iCone) = 0; % H1 hor. cells receive little S-cone input: cs x cs
	w = w ./ sum(w, 2); % sum of input weights is 1: cs x cs
	m.p.wConeHor = w; % cone to horizontal cell synaptic weights: cs x cs
	w = a'; w = w ./ sum(w, 2); % sum of input weights is 1: cs x cs
	m.p.wHorCone = w; % horizontal cell to cone synaptic weights: cs x cs

	% Calculate bipolar to off-ganglion cell attentuation
	locI = m.p.cone.loc; % bipolar cell locations: cs x 2
	locO = m.p.gangOff.loc; % ganglion cell locations: gs x 2
	radGang = m.p.radGang; % ganglion cell dendritic radius: 1 x 1
	w = calConv(locI, locO, radGang); % attenuation for off-cells: gFs x cs
	w = w ./ sum(w, 2); % sum of input weights is 1: gFs x bs
	m.p.wBiGangOff = w; % bipolar cell to off-gang. synaptic weights: gFs x cs
	
	% Calculate bipolar to on-ganglion cell attentuation	
	locO = m.p.gangOn.loc; % ganglion cell locations: gs x 2
	w = calConv(locI, locO, radGang); % attenuation for on-cells: gNs x cs
	w(:, iCone) = 0; % there are no on-S midget ganglion cells: gNs x bs
	w = w ./ sum(w, 2); % sum of input weights is 1: gNs x bs
	if radGang == 0 % private line
		i = isnan(w); % nan results from dividing 0 by 0
		w(i) = 0; % weights are 0
	end
	m.p.wBiGangOn = w; % bipolar cell to on-gang. synaptic weights: gNs x bs

	% Calculate geniculate to 4CBeta excitatory cell attentuation
	locI = [m.p.gangOff.loc; m.p.gangOn.loc]; % geniculate cell locations: is x 2
	locO = m.p.betaEx.loc; % 4CBeta excitatory locations: os x 2
	r = m.p.radBeta; % geniculocortical convergence radius (deg)
	a = calConv(locI, locO, r); % attenuation: os x is
	if radGang == 0 % private line
		a(:, iCone) = 0; % off-S geniculate cells do not project to 4CBeta
	end
	a = a ./ sum(a, 2); % sum of input attenuations is 1: os x is
	m.p.a_gen_betaEx = a; % geniculocortical attenuation: os x is

	% Calculate geniculate to 4CBeta inhibitory cell attentuation
	locO = m.p.betaIn.loc; % 4CBeta inhibitory locations: os x 2
	a = calConv(locI, locO, r); % attenuation: os x is
	if radGang == 0 % private line
		a(:, iCone) = 0; % off-S geniculate cells do not project to 4CBeta
	end
	a = a ./ sum(a, 2); % sum of input attenuations is 1: os x is
	m.p.a_gen_betaIn = a; % geniculocortical attenuation: os x is

	%	Increment off-centre geniculocortical gain to stop private line cancellation
	if m.p.varyGain % flag is on
		%	modC = normrnd(1, m.p.kDevGain, [1, size(locI, 1)]); % deviates: 1 x is
		modC = ones(1, size(locI, 1)); % modulation for all inputs: 1 x is
		modC(1: size(m.p.gangOff.loc, 1)) = m.p.kGainOff; % inc. for offs: 1 x is
		m.p.a_gen_betaEx = modC .* m.p.a_gen_betaEx; % multiply
		m.p.a_gen_betaIn = modC .* m.p.a_gen_betaIn; % multiply
	end		

	% Calculate 4CBeta inhibitory to excitatory attenuation
	locI = m.p.betaIn.loc; % 4CBeta inhibitory cell locations: is x 2
	locO = m.p.betaEx.loc; % 4CBeta excitatory cell locations: os x 2
	r = m.p.radInEx; % inhibitory-excitatory convergence radius (deg)
	a = calConv(locI, locO, r); % attenutation: os x is
	a = a ./ sum(a, 2); % sum of input weights is 1: os x is
	m.p.w_betaIn_betaEx = a; % inhibitory to excitatory attenuation: os x is

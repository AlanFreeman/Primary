%	Library of task functions

function h = runPrimTask % obtain handles of task functions in this file

	h = localfunctions; % handles of task functions in this file
		
function [d, m] = addW(d, m) % add synaptic weight to stage properties

	%	Calculate synaptic weights
	d.prop.weight = []; % allocate storage
	s = d.stage; % name of current stage
	i = ismember(s, {'genOff', 'genOn'}); % only eligible stages
	if i == 0, return, end % nothing to do
	m = m.setWeight(d, m); % calculate weights

	%	Store synaptic weights
	for post = ["betaEx", "betaIn"] % loop over postsynatic stages
		w = "w_gen_" + post; % name of current synaptic weight
		w = m.p.(w); % value of current synaptic weight: os x is
		loc = m.p.(post).loc; % locations of postsynaptic cells (deg): ls x 2
		i = knnsearch(loc, m.addW.loc); % index of specified postsynaptic cell
		w = w(i, :); % weight: all presyn. cells with spec. postsyn. cell: 1 x is
		i = size(m.p.gangOff.loc, 1); % number of off-centre geniculate cells
		if s == "genOff", w = w(1: i); % offs come first in weight array
		else, w = w(i + 1: end); % on-geniculates
		end
		d.prop.weight.(post) = w; % store
	end

function [d, m] = centre(d, m) % correct resp. phase for distance from origin

	% Find the model parameters to be varied
	[name, val, vals] = m.getVal(m, 'resp'); % model parameter names and values
	vs = prod(vals); % total number of values
	r = d.resp; % response: 1 x fs x ls x 1 x vs, vs can be multi-dimensional
	s = size(r); sC = [s(1: 4), vs]; % vs is 1-dimenional: 1 x fs x ls x 1 x vs
	r = reshape(r, sC); % reshape the response

	% Determine response phase shift for stimulus centred on cell
	for v = 1: vs % loop over variable values

		% Calculate distance of cell from origin
		m = m.setVal(m, name, val, v); % set values of model parameters
		m = m.samp(m); % recalculate temporal par. if temp. freq.	changes
		l = d.loc; % cell location (deg): 1 x ls x 2
		l = shiftdim(l, 1); % make it a matrix: ls x 2
		dir = pi * m.p.dir / 180; % grating direction (radians)
		freqS = 2 * pi * m.p.freqS; % grating spatial frequency (radians/deg)
		u = [cos(dir); sin(dir)]; % unit vector in grating direction: 2 x 1
		p = l * u; % distance of cell from origin in grating direction (deg): ls x 1

		% Correct phase
		p = freqS .* p; % phase shift from origin (radians): ls x 1
		p = exp(1i * p); % phase correction (complex): ls x 1
		p = shiftdim(p, -2); % prepare to correct phase: 1 x 1 x ls
		r(:, :, :, :, v) = p .* r(:, :, :, :, v); % correct phase: 1 x fs x ls x vs

	end

	% Store
	r = reshape(r, s); % restore response shape: vs can be multi-dimensional
	d.resp = r; % response with phase corrected: 1 x fs x ls x vs

function [d, m] = comp(d, m) % calc. cone components in centre and sur. resp.

	r = d.resp; % fundamental for all cone types (mV): 1 x 1 x ls x zs x cs
	r = abs(r); % response magnitude (mV): 1 x 1 x ls x zs x cs
	r = r ./ sum(r, 5); % normalise by sum over cone types: 1 x 1 x ls x zs x cs
	d.resp = r; % normalised cone components: 1 x 1 x ls x zs x cs
	d.Properties.VariableDescriptions{'resp'} = 'Cone ratio'; % describe

function [d, m] = corCourse(d, m) % calc. cross-correlation between time courses

%	Assumption: time courses are in the frequency domain

	%	Prepare for cross-correlation
	locI = d.resp(1).loc; % location (deg): ls x 2
	rI = d.resp(1).resp; % input time course, frequency domain (Hz): ts x is
	rI(1, :) = 0; % subtract mean: ts x is
	rI = permute(rI, [1, 3, 2]); % ts x 1 x is
	locO = d.resp(2).loc; % location (deg): ls x 2
	rO = d.resp(2).resp; % output time course, freq. domain (Hz): ts x os
	rO(1, :) = 0; % subtract mean: ts x os

	%	Cross-correlate
	c = conj(rI) .* rO; % transform of cross-correlation (Hz^2): ts x os x is
	c = ifft(c); % cross-correlation, time domain (Hz^2): ts x os x is
	c = c ./ max(c); % normalised cross-correlation: ts x os x is

	%	Store
	d.locI = shiftdim(locI, -1); d.locO = shiftdim(locO, -1); % store: 1 x ls x 2
	d.resp = shiftdim(c, -1); % store: 1 x ts x os x is
	d.Properties.CustomProperties.RespDim = {'', 'time', 'locO', 'locI'}; % update

function [d, m] = crop(d, m) % crop locations to smaller field

	% Prepare response
	r = d.(m.z); % response: ps x ls x qs; ps and qs can be multi-dimensional
	loc = shiftdim(d.loc, 1); % locations (deg): ls x 2
	if isfield(m.crop, 'radius') % radius is specified
		rad = m.crop.radius; % radius (deg)
	else % default radius
		rad = m.p.radSur; % surround radius (deg)
	end
	rad = .5 * m.p.wid - rad; % half-width of cropped field (deg)
	if isfield(m.crop, 'dim') % set location dimension
		dim = m.crop.dim; % user specified
	else % default
		dim = d.Properties.CustomProperties.RespDim; % response dimensions
		[~, dim] = ismember('loc', dim); % index of location in response
	end

	%	Crop and store
	[r, loc] = m.doCrop(r, loc, rad, dim);
	d.(m.z) = r; % store response: ps x ls x qs, ps and qs can be multi-dim.
	d.loc = shiftdim(loc, -1); % locations (deg): ls x 2

function [d, m] = dens(d, m) % calculate cone density

	[densC, x] = prim.densFunc(m); % cone density cross-cor. with conv. function
	d.x = x; d.y = x; % x and y locations (deg)
	d.dens = shiftdim(densC, -1); % store
	d.Properties.VariableDescriptions{'dens'} = 'Cone density (cones / v.f.)';

function [d, m] = dep(d, m) % select one of multiple values in a dependent var.

	r = d.(m.z); % dependent variable
	[r, dim] = prim.readDep(r, m.dep.name); % select single value from dep't var.
	s = size(r); % size of (new) dependent variable
	if s(1) ~= 1 % cannot store variable with more than one row
		r = shiftdim(r, -1); % set to one row
	end
	d.(m.z) = r; % store
	dim = [{''}, dim]; % dimensions
	d.Properties.CustomProperties.RespDim = dim; % reset dimensions
	
function [d, m] = desc(d, m) % describe the variables in the data table

	desc = d.Properties.VariableDescriptions; % variable descriptions
	undesc = find(strcmp('', desc)); % indices of undescribed variables
	for i = undesc % loop over undescribed variables
		name = d.Properties.VariableNames{i}; % name of current variable
		switch name % describe this variable
			case 'cont', desc{i} = 'Contrast';	
			case 'dir', desc{i} = 'Stimulus direction (deg)';
			case 'freqS', desc{i} = 'Spatial frequency (cycles/deg)';
			case 'freqSPref', desc{i} = 'Preferred spatial frequency (cycles/deg)';
			case 'freqT', desc{i} = 'Temporal frequency (Hz)';
			case 'loc', desc{i} = 'Visual field location (deg)';
			case 'resp', desc{i} = 'Generator potential (mV)';
			case 'stage', desc{i} = 'Processing stage';
			case 'time', desc{i} = 'Time (s)';
			case 'x', desc{i} = 'Horizontal location (deg)';
			case 'xs', desc{i} = 'Number of locations';
			case 'y', desc{i} = 'Vertical location (deg)';
		end
	end
	d.Properties.VariableDescriptions = desc; % store

function [d, m] = doMax(d, m) % find maximum response and tuning curves

	%	Initialise
	z = d.(m.z); % response (mV or Hz): 1 x 1 x ls x 1 x vs, vs may be multi-dim.
	ls = size(z, 3); % number of locations
	par = prim.getPar(d, m); % stimulus parameters and their values
	name = [par.name]; % parameter names: 1 x vs, string
	vs = length(name); % number of stimulus parameters
	z = reshape(z, [ls, par.vals]); % remove single-element dim.: ls x vs

	% Calculate maximum, and organise indices
	[zMax, iMax] = max(z, [], 1 + (1: vs), 'linear'); % maxima and their indices
	sub = cell(1, 1 + vs); % allocate storage
	[sub{:}] = ind2sub(size(z), iMax); % indices of pref. stimulus values
	sub = sub(1 + (1: vs)); % discard first element (which is 1: ls): 1 x vs
	sub = cell2mat(sub); % convert indices to matrix: ls x vs

	%	Obtain tuning curves of stimulus variables
	zTun = struct('val', []); % allocate storage
	for i = 1: vs % loop over stimulus variables
		ord = 1: vs; ord(i) = []; % permutation order
		zP = permute(z, 1 + [0, i, ord]); % permute z to place variable i second
		subC = sub; subC(:, i) = []; % subscripts other than variable i
		for j = 1: ls % loop over locations
			k = subC(j, :); % indices for this location
			k = num2cell(k); % convert to cell
			zTun(i).val(j, :) = zP(j, :, k{:}); % tuning curve: ls x ts
		end
	end

	%	Store
	for i = 1: vs % loop over stimulus variables
		nC = name{i}; % name of current variable
		resp.(nC).pref.val = d.(nC)(sub(:, i))'; % stimulus preference: ls x 1
		resp.(nC).pref.dim = {'loc'}; % dimensions
		resp.(nC).tun.val = zTun(i).val; % stimulus tuning: ls x ts
		resp.(nC).tun.dim = {'loc', nC}; % dimensions
	end
	d.resp = resp; % preferences and tuning
	d.resp.max.val = zMax(:); % maximised response: ls x 1
	d.resp.max.dim = {'loc'}; % dimensions
	d.resp.full.val = z; % full response: ls x vs
	d.resp.full.dim = ['loc', par.name]; % dimensions

function [d, m] = doMod(d, m) % calculate synaptic modulation factors

	% Initialise development data
	switch m.mod.init % initialise or continue
		case 'cont' % continue from existing file
			phase = 2; % development phase
			syn = d.syn(end, :); % plastic synaptic sites at end of previous phase
			d(end + 1, :) = d(end, :); % 0th cycle of phase 2 = final cycle of phase 1
			d.file(end) = m.mod.serial; % set file number
			d.phase(end) = 2; d.cycle(end) = 0; % set phase, cycle
		case 'init' % initialise file

			%	Variables to store
			phase = 1; % development phase
			locGen = [m.p.gangOff.loc; m.p.gangOn.loc]; % geniculate cell loc. (deg)
			locGen = shiftdim(locGen, -1); % for storage: 1 x ls x 2
			locBetaEx = shiftdim(m.p.betaEx.loc, -1); % ex. loc. (deg): 1 x os x 2
			locBetaIn = shiftdim(m.p.betaIn.loc, -1); % in. loc. (deg): 1 x os x 2	
			syn = m.p.syn; % plastic synaptic sites

			% Create file
			d = table(m.mod.serial, m.p.ecc, m.p.wid, phase, 0, locGen, ...
				locBetaIn, locBetaEx, syn, 'variableNames', {'file', 'ecc', 'width', ...
				'phase', 'cycle', 'locGen', 'locBetaIn', 'locBetaEx', 'syn'});
			d.Properties.Description = m.p.project; % research project
			d.Properties.VariableDescriptions{1} = ''; % set all descriptions to empty

	end
	cycles = max(m.p.cycle); % number of development cycles
	kInc = m.p.kIncMod / cycles; % maximum change of modulation factor per cycle
	syns = size(syn, 2) - 1; % number of excitatory sites

	% Loop over development cycles and update weights on each cycle
	[name, val, vals] = m.getVal(m, 'mod'); % stimulus names and values
	vs = prod(vals); % total number of stimuli
	dS = cell(length(m.p.cycle), 1); % allocate storage for data files
	for cycle = 1: cycles % loop over cycles

		%	Loop over synaptic sites and install synaptic weight for each site
		for i = 1: syns % loop over synaptic sites
			att = syn(i).att; % name of distance-based attenuation
			w = syn(i).mod .* m.p.(att); % synaptic weights: os x is
			w = w ./ sum(w, 2); % normalise to unity sum: os x is
			m.p.(syn(i).weight) = w; % store: os x is
			sizeC = size(w); % size of synaptic weight matrix
			syn(i).inc = zeros([sizeC, vs]); % modulation increment: os x is x vs
		end
		m.p.kInEx = syn(end).mod; % inhibitory-excitatory gain: 1 x 1

		% Loop over stimuli and calculate cross-correlations
		ls = size(m.p.betaEx.loc, 1); % number of excitatory cells
		r = zeros(ls, vs); % maximum excitatory response, to calculate inhib. plast.
		for v = 1: vs % loop over stimuli

			% Calculate input and output time courses
			m = m.setVal(m, name, val, v); % set values of model parameters
			pAll = m.solveF(m); % frequency-domain response (mV): ls x ts

			for i = 1: syns % calculate modulation factor increments and store them

				%	Extract input and output stages
				s = syn(i).stage; % names of input and output stages
				if s(1) == "gen", s = ["genOff", "genOn", s(2)]; end % expand 'gen'
				j = contains({pAll.type}, s); % required stages
				p = vertcat(pAll(j).resp); % response (mV): ls x ts
				p = p .'; % put time dimension first: ts x ls
	
				% Calculate cross-correlation between input and output for each synapse
				p = m.rect(p); % potential proportional to impulse rate (mV): ts x ls
				p(1, :) = 0; % subtract mean (mV): ts x ls
				is = size(syn(i).inc, 2); % number of inputs
				pS = p(:, 1: is); % presynaptic input (mV): ts x is
				pC = p(:, is + 1: end); % postsynaptic input (mV): ts x os
				syn(i).inc(:, :, v) = m.cor(m, pS, pC); % cor'n (mV^2): os x is x vs

			end

			%	Store excitatory response for calculation of inhibitory plasticity
			i = contains({pAll.type}, "betaEx"); % index of excitatory respons
			p = pAll(i).resp; % excitatory potential (mV): ls x ts
			p = ifft(p, [], 2); % convert to time domain (mV): ls x ts
			p = real(p); % remove any imaginary components: ls x ts
			p = m.p.kRect * max(0, p); % impulse rate (Hz): ls x ts
			r(:, v) = max(p, [], 2); % maximum impulse rate (Hz): ls x vs

		end

		% Adjust increments and add to modulation
		%	rad = .5 * m.p.wid - m.p.radBetaSur; % v.f. radius to remove edge effects
		for i = 1: syns % loop over synaptic sites
			inc = syn(i).inc; % increment for this site: os x is x vs
			[~, j] = max(abs(inc), [], 3, 'linear'); % index of extreme cor.: os x is
			inc = inc(j); % remove suboptimal correlations: os x is
			inc = real(inc); % remove any imaginary components: os x is
			%{
			s = syn(i).stage; % stages involved in current synapse
			locBeta = m.p.(s(2)).loc; % locations of postsynaptic cells: ls x 2
			incC = m.doCrop(inc, locBeta, rad, 1); % exclude edge effects
			inc = inc / (2 * std(incC(:))); % normalise by near maximum: os x is
			%}
			inc = inc / (2 * std(inc(:))); % normalise by near maximum: os x is
			inc = kInc * inc; % set maximum decrement and increment: os x is
			modC = syn(i).mod + inc; % update modulation factor: os x is
			syn(i).mod = min(2, max(0, modC)); % limit modulation factor: os x is
		end

		%	Reset inhibitory-excitatory gain
		r = max(r, [], 2); % maximum response across stimuli (Hz): ls x 1
		%	r = m.doCrop(r, m.p.betaEx.loc, rad, 1); % exclude edge effects: ls x 1
		r = max(r); % maximum response across excitatory cells (Hz): 1 x 1
		syn(end).mod = 1 + m.p.kIncIn * r; % gain is proportional to max. response

		% Store
		[t, i] = ismember(cycle, m.p.cycle); % store this cycle?
		if t % yes
			dC = d(end, :); % single-row table in which to store current cycle
			dC.file = m.mod.serial; % modulation file serial number
			dC.phase = phase; dC.cycle = cycle; % update
			syn = rmfield(syn, 'inc'); % interim result not required for storage
			dC.syn = syn; % update modulations
			if isfield(m.mod, 'save') % save interim file
				filename = "Modulation cycle " + int2str(cycle); % file name
				save(filename, 'dC'); % save
			end
			dS{i} = dC; % store for concatenation
		end

	end
	d = vertcat(d, dS{:}); % combine data tables

function [d, m] = doPca(d, m) % principal components analysis

	out = prim.pcaFunc(d.resp, m); % perform PCA
	name = string(fieldnames(out))'; % field names of out
	for n = name % loop over field names
		d.resp.(n).val = out.(n); % store in table
	end

function [d, m] = field(d, m) % calculate receptive field by cross correlation

% Method: deliver a range of pulsed stimuli, weight each by response, and add

	%	Set the values of the stimulus variables to be varied
	name = d.Properties.CustomProperties.RespDim(5: end); % names of variables
	names = length(name); % number of variables
	val = cell(1, names); vals = zeros(1, names); % values, number of values
	for i = 1: names % loop over variables
		val{i} = d.(name{i})'; % assign values
		vals(i) = length(val{i}); % number of values
	end

	% Calculate the stimulus pixel locations	
	ss = prod(vals); % total number of stimuli
	switch m.p.stimS % grating or image
		case 'grating'
			loc = m.listLoc(m); % locations (x, y) (deg):	ps x 2
			ps = size(loc, 1); % number of stimulus pixels
		case 'image'
			[xs, ys, ~] = size(val); % image size
			ps = xs * ys; % number of stimulus pixels
	end
	
	% Obtain the responses
	r = d.resp; % impulse rate (Hz): 1 x ts x ls x 1 x ss, ss can be multi-dim.
	r = reshape(r, [], ss); % prepare for cross correlation: ls x ss
	rSum = sum(r, 2); % sum of responses over all stimuli: ls x 1

	%	Cross correlate
	s = zeros(ss, ps * 3); % stimuli: ss x pcs, pcs = ps * cs
	for i = 1: ss % loop over stimuli
		m = m.setVal(m, name, val, i); % set values of model parameters
		switch m.p.stimS % grating or image
			case 'grating'
				if ismember('cont', d.Properties.VariableNames) % cont. is a variable
					m.p.cont = d.cont; % set it
				end
				stim = m.calStim(0, loc, m); % stimulus (contrast-units): 1 x ps x cs
			case 'image'
				stim = m.p.image; % current map: xs x ys x cs
				stim = m.rgb2cont(m, stim); % convert to contrast units
				stim = reshape(stim, 1, ps, []); % stimulus (cont.-units): 1 x ps x cs
		end
		s(i, :) = stim(1, :); % prepare for multiplication: ss x pcs
	end
	f = r * s; % receptive field (Hz x c.-u.): (ls x ss) * (ss x pcs) = ls x pcs
	f = reshape(f, [], ps, 3); % reshape: ls x ps x cs
	%	f = f / (ss * m.p.contMag); % normalise (Hz): ls x ps x cs
	f = f ./ (rSum * m.p.contMag); % normalise: ls x ps x cs

	%	Store
	d.locS = shiftdim(loc, -1); % stimulus locations (deg): 1 x ps x 2
	d.colour = {'L', 'M', 'S'}; % colour space
	d.resp = shiftdim(f, -1); % receptive field (Hz): 1 x ls x ps x cs
	d.Properties.CustomProperties.RespDim = ...
		{'', 'loc', 'locS', 'colour'}; % response dimensions

function [d, m] = fig(d, m) % calculate functions for figure plot

	n = fieldnames(m.fig)'; % cell array of names including function types: 1 x ts
	for fun = string(n) % loop over function types
		switch fun % choose function type
			case 'conv' % convergence functions
	
				% Calculate convergence functions
				xs = 101;
				if isfield(m.fig, 'xLim') % x limit (deg)
					xLim = m.fig.xLim; % specified by user
				else, xLim = .1; % default
				end
				x = xLim * linspace(-1, 1, xs); % visual field locations (deg): 1 x xs
				rad = [m.p.radOpt; m.p.radHor; m.p.radGang; m.p.radBeta];
					% radii (deg): rs x 1
				f = exp(- x .^ 2 ./ rad .^ 2); % convergence function: rs x xs
			
				% Create a table for each radius
				source = ["opt", "hor", "gang", "beta"]; ss = length(source);
				dC = cell(1, ss); % tables
				for i = 1: ss % loop over sources
					dC{i} = table(source(i), x, rad(i), f(i, :), ...
						'variableNames', {'source', 'x', 'rad' 'f'}); % make table
				end
			
				% Concatenate tables
				d = vertcat(dC{:}); % concatenate
				d.source = categorical(d.source); % for listing
				d.Properties.VariableDescriptions{'source'} = 'Function source';
				d.Properties.VariableDescriptions{'x'} = 'Visual field location (deg)';
				d.Properties.VariableDescriptions{'rad'} = 'Function radius (deg)';
				d.Properties.VariableDescriptions{'f'} = 'Convergence function';
	
			case 'course' % stimulus time course
	
				s = d.stim; % stimulus: 1 x ts x xs x ys
				t = d.t; % times (s): 1 x ts
				i = t <= 1 / m.p.freqT; % reduce to first cycle (s): 1 x ts, new ts
				x = d.x; y = d.y; % locations (deg): 1 x xs
				xs = length(x); % number of locations
				[x, y] = ndgrid(x, y); % location grids (deg): xs x ys
				loc = [x(:), y(:)]; % locations (deg): xs * ys x 2
				j = knnsearch(loc, [0, 0]); % location nearest to middle of visual field
				[j, k] = ind2sub([xs, xs], j); % middle location (deg): 1 x 2
				s = s(1, i, j, k); % time course of stimulus at middle: 1 x ts
				d.t = t(i); % first cycle only (s): 1 x ts
				d.course = s; % store: 1 x ts
	
			case 'mech' %  % circles with radii equal to centre, surround
				loc = m.fig.mech; l = [loc; loc]; % mechanism locations
				r = [m.p.radCen; m.p.radSur]; % mechanism radii
				viscircles(l, r); % draw mechanisms
		end
	end

function [d, m] = fitGauss(d, m) % fit Gaussian to spatial frequency response

	b(1) = d.resp(1); % peak response
	rCen = sqrt(m.p.radOpt ^ 2 + m.p.radGang ^ 2); % estimated centre rad. (deg)
	rSur = sqrt(2 * m.p.radHor ^ 2); % estimated surround radius (deg)
	switch m.fitGauss.comp % centre or surround?
		case 'cen', r = rCen; % centre component
		case 'sur', r = rSur; % surround component
	end
	switch m.fitGauss.coef % number of regression coefficients
		case 1 % peak response only
			fun = @(b, f) b * exp(-.25 * (r * 2 * pi * f) .^ 2); % Gaussian
		case 2 % peak response and radius
			b(2) = r; % radius (deg)
			fun = @(b, f) b(1) * exp(-.25 * (b(2) * 2 * pi * f) .^ 2); % Gaussian
		case 4 % both centre and surround components
			b(2) = rCen; b(3) = .5 * b(1); b(4) = rSur; % estimate coefficients
			fun = @(b, f) b(1) * exp(-.25 * (b(2) * 2 * pi * f) .^ 2) - ...
				b(3) * exp(-.25 * (b(4) * 2 * pi * f) .^ 2); % difference of Gaussians
	end
	model = fitnlm(d, fun, b, 'predictorVars', {'freqS'}, 'responseVar', 'resp');
		% fit data
	m.model{m.group} = model;
	
function [d, m] = fund(d, m) % calculate fundamental Fourier response

% Assumptions:
%		frequency is the second dimension of d.resp

	% Calculate fundamental
	r = d.resp; % response: rs x fs x vs where vs can be multi-dimensional
	s = size(r); % response size
	rs = s(1); fs = s(2); % number of rows, frequencies
	r = r(:, 2, :); % keep only fundamental component: rs x 1 x vs, vs single-dim.
	r = reshape(r, [rs, 1, s(3: end)]); % restore shape: rs x 1 x vs, vs multi.
	r = (2 / fs) * r; % convert from fft units to temporal units: rs x 1 x vs
	
	% Store
	if isfield(m.fund, 'prop') % response property: amplitude, complex or phase
		p = m.fund.prop; % user-specified
	else
		p = 'amp'; % default: amplitude
	end
	switch p % amplitude, complex or phase?
		case 'amp' % amplitude
			d.resp = abs(r); % store amplitude: rs x fs x vs
			%	desc = 'Fundamental amplitude'; % response description
		case 'complex'
			d.resp = r; % store: rs x fs x vs
			%	desc = 'Fundamental component';
		case 'phase'
			d.resp = (180 / pi) * angle(r); % store phase (deg): rs x fs x vs
			d.Properties.VariableDescriptions{'resp'} = 'Fundamental phase (deg)';
	end
	d.freq = d.freq(:, 2); % update frequency

function [d, m] = hist(d, m) % calculate histogram

	% Initialise
	if isfield(m.hist, 'edges') % set edges or number of bins
		arg = m.hist.edges; % user-specified
	else
		arg = 10; % number of bins
	end
	
	% Calculate histogram
	x = d.(m.x); % data values
	[n, e] = histcounts(x, arg); % obtain histogram counts and edges
	w = e(2) - e(1); % bin width
	c = .5 * w + e(1: end - 1); % bin centres
	d = d(1, :); % keep only the first row
	d.(m.x) = c; % bin centres
	d.edge = e; % bin edges
	d.(m.z) = n; % counts
	d.Properties.VariableDescriptions{'edge'} = 'Edge';
	d.Properties.VariableDescriptions{m.z} = 'Count';

function [d, m] = index(d, m) % calculate response indices

%	Input:
%	*	some or all stimulus types: M-specific, equiluminant, L-specific, achromatic
%	*	outputs from task "max"
%	Output:
%	*	d.resp = structure with one field for each index, e.g. coi, soi

	%	Calculate indices
	if isfield(m.index, "name") % is an index required?
		n = string(m.index.name); % indices to calculate: 1 x is, string
		index = prim.indexFunc(d, m, n); % calculate indices
	else, return % nothing to do
	end

	%	Store
	d = d(1, :); % keep just the first row
	d.resp = struct; % initialise as structure
	for nC = n % loop over indices
		d.resp.(nC).val = index.(nC); % chromatic opponency index: ls x vs
		d.resp.(nC).dim = {'loc'}; % response dimensions
	end

function [d, m] = interp(d, m) % fit tuning curve with an interpolated model

	% Calculate x values for interpolation
	x = d.(m.x); % tuning variable (deg or cycles/deg): 1 x xs
	z = d.(m.z); % response: structure or 1 x ls x xs
	if isstruct(z), z = z.val; else, z = shiftdim(z, 1); end % response: ls x xs
	z = z'; % transpose for interp1: xs x ls
	if m.x == "dir" % cyclic stimulus?
		cyclic = 1; % yes
	else, cyclic = 0; % no
	end
	if cyclic % stimulus is cyclic: add starting point at end, for interp.
		x(end + 1) = x(1) + 360; % close open-ended interval (deg): xs -> xs + 1
		z(end + 1, :) = z(1, :); % cyclic z: xs x ls
	end
	if isfield(m.interp, 'xs') % set number of interpolation points
		xs = m.interp.xs; % user-defined
	else, xs = 101; % default
	end
	xI = linspace(x(1), x(end), xs); % interpolation points: 1 x xs

	% Interpolate, find maximum and preferred direction
	if isfield(m.interp, 'method'), method = m.interp.method; % set method
	else, method = 'spline'; end % default
	z = interp1(x, z, xI', method); % interpolate: xs x ls
	[zMax, i] = max(z); % maximum response, and its index: 1 x ls
	xPref = xI(i); % preferred value of tuning curve: 1 x ls

	%	Calculate tuning curve bandwidth
	zCrit = zMax / sqrt(2); % response criterion used by Ringach (02): 1 x ls
	w = m.bandwidth(xI, z, zCrit, cyclic); % left- and right-bandwidths: 2 x ls
	w = mean(w); % mean bandwidth: 1 x ls

	%	Store
	%	d.(m.x) = xI; % tuning variable: 1 x xs
	zC = permute(z, [3, 2, 1]); % response: 1 x ls x xs
	if isfield(m.interp, 'out') % choose output variable
		z = m.interp.out; % user-specified name
	else, z = m.z; % default: same as input name
	end
	switch z
		case 'band' % bandwidth
			d.band = w; % bandwidth: 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'dirPref' % preferred direction
			d.dirPref = xPref; % preferred direction (deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'dirTun' % direction tuning
			d.(m.x) = xI; % tuning variable: 1 x xs
			d.dirTun = zC; % direction tuning: 1 x ls x ds
			dim = {'', 'loc', 'dir'}; % update output dimensions
		case 'freqSPref' % preferred spatial frequency
			d.freqSPref = xPref; % preferred spatial freq. (cycles/deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'freqSTun' % spational frequency tuning
			d.(m.x) = xI; % tuning variable: 1 x xs
			d.freqSTun = zC; % store spatial frequency resp.: 1 x ls x fs
			dim = {'', 'loc', 'freqS'}; % update output dimensions			
	end
	d.Properties.CustomProperties.RespDim = dim; % update output dimensions

function [d, m] = loc(d, m) % find location of cell plotted on existing axes

%	Input: scatter plot with cell locations in UserData, and cell datatipped

	loc = prim.see(m, locate = 'scatter'); % cell location (deg): 1 x 2
	fprintf('Location of datatipped cell (deg): %.4g, %.4g\n', loc(1), loc(2));

function [d, m] = mark(d, m) % mark points on existing axes

%	Inputs:
%		m.mark.loc = locations of cells to mark [x1, y1; x2, y2; ...] (deg): ms x 2
%		locAll = name of variable containing locations of all cells (deg): ls x 2
%	Assumptions:
%		there is one point for each cell

	%	Find indices of points to mark
	if isfield(m.mark, "loc") % any points to mark?
		loc = m.mark.loc; % cell locations [x1, y1; x2, y2; ...] (deg): ms x 2
	else
		return % nothing to do
	end
	if isfield(m.mark, "locAll") % locations of all cells (deg): 1 x ls x 2
		locAll = d.(m.mark.locAll); % set by user
	else
		locAll = d.loc; % default
	end
	if ndims(locAll) == 3 % location variable has form 1 x ls x 2
		locAll = shiftdim(locAll, 1); % remove leading dimension
	end
	i = knnsearch(locAll, loc); % index of locations to mark

	%	Mark
	a = gca; % current axes
	a.NextPlot = "add"; % add to existing axes, don't replace
	o = a.Children; % object in which to mark, e.g. scatter plot
	x = o.XData(i); % x values of point to mark
	y = o.YData(i); % y values of point to mark
	plot(a, x, y, "ok", LineWidth = 2); % mark points with circle

function [d, m] = move(d, m) % show a movie

	%	Initialise
	im = shiftdim(d.image, 1); % image: ys x xs x cs x fs
	fs = size(im, 4); % number of frames
	x = d.x; y = d.y; % x and y limits (deg): 1 x 2
	dur = d.time(2) - d.time(1); % frame duration (s)
	a = gca; % current axes
	a.NextPlot = 'add'; % retain axis settings during movie

	%	Show movie
	for i = 1: fs % loop over frames
		imC = im(:, :, :, i); % current image: ys x xs x cs
		image(a, x, y, imC); % draw image
		pause(dur); % pause between frames (s)
	end

function [d, m] = op(d, m) % calculate opponency index

	%	Calculate
	fL = shiftdim(d.resp(2, :, :, 1), 1); % L-specific receptive field: ls x ps
	fM = shiftdim(d.resp(1, :, :, 2), 1); % M-specific receptive field: ls x ps
	ind = m.oi(fL, fM); % calculate opponency index

	%	Store
	d = d(1, :); % keep first row only
	d.oi = ind'; % opponency index, [0, 1]: 1 x ls
	d.Properties.CustomProperties.RespDim = {'', 'loc'}; % response dimensions

function [d, m] = peak(d, m) % find time course maximum

% Assumptions:
%		time is the second dimension of d.resp
%		d contains the variable 'time'

	r = d.resp; % response: rs x ts x vs where rs = rows, vs can be multi-dim.
	[d.resp, i] = max(r, [], 2); % peak response: rs x 1 x vs
	t = d.time(i(:)); % peak times (s): 1 x vs
	d.time = reshape(t, size(d.resp)); % update time: rs x 1 x vs

function [d, m] = predOr(d, m) % calc. pref. orient. from g.c. coverage map

	% Prepare coverage map
	xs = d.xs; % number of x locations
	if isfield(m.orient, 'debug') % debug using a cosine grating
		dir = (pi / 180) * m.p.dir; % motion dir. (radians ACW from rightward)
		fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
		pS = (pi / 180) * m.p.phaseS; % spatial phase (radians)
		x = .5 * m.p.wid * linspace(-1, 1, xs); % x values (deg): 1 x xs
		[x, y] = ndgrid(x); % grid locations (deg): xs x ys
		u = cos(dir) * x + sin(dir) * y; % dist. relative to centre (deg): xs x ys
		z = cos(fS * u - pS); % stimulus: xs x ys
	else % use coverage map
		z = d.(m.z); % coverage map: 1 x ls
		z = reshape(z, [xs, xs]); % prepare for imgaborfilt: xs x ys
	end
	z = z'; % make it an image: ys x xs
	z = flipud(z); % y increasing downwards: ys x xs

	% Make Gabor filter bank
	w = m.orient.freqS; % spatial frequency (cycles/deg): 1 x fs
	w = (xs / m.p.wid) ./ w; % wavelengths (pixels/cycle): 1 x fs
	o = m.orient.dir; % imgaborfilt-orientation = user-direction (deg): 1 x os
	g = gabor(w, o); % filter bank: 1 x gs, gs = fs * os
	w = [g.Wavelength]; % filter bank wavelength: 1 x gs
	f = (xs / m.p.wid) ./ w; % filter bank spatial frequency (cycles/deg): 1 x gs
	o = [g.Orientation]; % filter bank orientation (deg): 1 x gs
	o = mod(o + 90, 180) - 90; % orientation (deg), [-90, 90): 1 x gs

	% Apply Gabor filters and find filters with maximum magnitude
	z = imgaborfilt(z, g); % mag. of Gabor-filtered coverage map: ys x xs x gs
	[~, i] = max(z, [], 3); % maximum Gabor magnitude over filters: ys x xs
	i = i(:); % make index a column vector: ys * xs x 1
	f = f(i); f = reshape(f, [xs, xs]); % pref. spat. freq. (cycles/deg): ys x xs
	o = o(i); o = reshape(o, [xs, xs]); % preferred orientation: ys x xs

	% Store
	d.x = .5 * m.p.wid * [-1, 1]; d.y = d.x; % location limits (deg)
	d.freqS = m.orient.freqS; % spatial frequency (cycles/deg): 1 x fs
	d.orient = m.orient.dir; % direction (deg): 1 x os
	d.freqSPref = shiftdim(f, -1); % preferred spatial frequency: 1 x ys x xs
	d.orPref = shiftdim(o, -1); % pref. orientation (deg), [-90, 90): 1 x ys x xs
	d.Properties.CustomProperties.RespDim = {'', 'y', 'x'}; % update

function [d, m] = prep(d, m) % prepare for plot

	%	Prepare
	if isfield(m.prep, 'prep') % is m.prep.prep defined?
		n = m.prep.prep; % preparations, cell array: 1 x ps
	else, return % nothing to do
	end
	for p = string(n) % loop over preparations
		switch p % choose preparation
			case 'cont2rgb' % convert contrast to RGB *** remove ***
				r = d.(m.z); % weighted stimulus
				%	r = r / max(r(:)); % weighted stimulus contrast, range is [-1, 1]
				d.(m.z) = .5 * (1 + r); % shift range to [0, 1]
			case 'image' % prepare for plotting with image or imagesc

			%	Assumptions:
			%		x values are located in d.(m.x)
			%		z values are located on a square array, d.(m.z)
			%		one of the RespDim dimensions is m.x

				%	Obtain pixel locations and response dimensions
				x = shiftdim(d.(m.x), 1); % pixel locations (deg): ls x 2
				xs = sqrt(size(x, 1)); % number of x values
				dim = d.Properties.CustomProperties.RespDim; % response dimensions
				i = contains(dim, m.x); % true for location dimension
				j = find(i); % index of location dimension
				iPre = 1: j - 1; iPost = j + 1: length(dim); % dim. before & after loc.
		
				% Reshape response
				z = d.(m.z); % image variable
				s = size(z); % image size
				z = reshape(z, [s(iPre), xs, xs, s(iPost)]); % map: ps x xs x ys x qs
				z = permute(z, [iPre, j + 1, j, iPost + 1]); % image: ps x ys x xs x qs

				%	Store
				%	xLim = [min(x(:, 1)), max(x(:, 1))]; % lower, upper bounds of x (deg)
				%	yLim = [min(x(:, 2)), max(x(:, 2))]; % lower, upper bounds of y (deg)
				%	d.x = xLim; d.y = yLim; % limits of x, y (deg): 1 x 2
				d.x = linspace(min(x(:, 1)), max(x(:, 1)), xs); % x (deg): 1 x xs
				d.y = d.x; % y (deg); contour plot requires all values: 1 x ys
				d.(m.z) = z; % store
				dim = [dim(iPre), 'x', 'y', dim(iPost)]; % update dimensions
				d.Properties.CustomProperties.RespDim = dim; % store

			case 'mark' % select cells from population
				loc = shiftdim(d.loc, 1); % cortical locations (deg): ls x 2
				i = knnsearch(loc, [m.p.c.b; m.p.c.c]); % indices of cells to select
				d.resp.max = d.resp.max(i); % maxima of selected cells (Hz): ls x 2
				d.resp.osi = d.resp.osi(i); % osi of selected cells: ls x 2
			case {'contour', 'mono', 'true'} % prepare for plotting receptive field
				switch p % monochromatic or true colour?
					case {'contour', 'mono'} % monochromatic
						if d.cont == m.p.contMag * [1, 1, 1] % achromatic stimulus
							d.resp = d.resp(:, :, :, :, 1); % all channels equal
						elseif d.cont == m.p.contMag * [1, -1, 0] % equiluminant stimulus
							d.resp = d.resp(:, :, :, :, 1); % rM = -rL, rS = 0
						elseif d.cont == m.p.contMag * [1, 0, 0] % L-specific stimulus
							d.resp = d.resp(:, :, :, :, 1); % L channel only
						else % M-specific stimulus
							d.resp = d.resp(:, :, :, :, 2); % M channel only
						end
					case 'true' % true colour
						r = d.(m.z); % response
						if d.cont == m.p.contMag * [1, -1, 0] % equiluminant stimulus
							d.(m.z) = .5 * (1 + r); % shift range to [0, 1]
						else % all other stimuli
							d = []; % delete
						end
				end
			case 'norm' % normalise by maximum
				v = d.(m.z); % values of variable: 1 x vs
				d.total = sum(v); % total number of values: 1 x 1
				maxC = max(v, [], 2); % 1 x 1
				d.(m.z) = v ./ maxC; % normalised values: 1 x vs
				d.Properties.VariableDescriptions{'total'} = 'Total count';
				d.Properties.VariableDescriptions{m.z} = 'Normalised count';
			case 'orient' % convert direction to orientation
				z = d.(m.z); % direction (deg), range = [-180, 180)
				z = z + 90; % range = [-90, 270)
				z = mod(z, 180); % range = [0, 180)
				d.(m.z) = z - 90; % orientation (deg), range = [-90, 90)
			case 'real' % make resp real *** check ***
				d.resp = real(d.resp); % make it real
			case 'resize' % resize image to, for example, smooth it *** fix ***

				% Obtain inputs
				if isempty(m.prep.resize) % scale is specified by user
					s = 4; % default: image width and height will be multiplied by s
				else, s = m.prep.resize; % user-specified
				end
			
				% Scale
				z = d.(m.z); % image: 1 x xs x ys x ps, ps can be multi-dimensional
				z = shiftdim(z, 1); % shift image to first two dimensions: xs x ys x ps
				xs = size(z, 1); % number of x values
				xs = s * (xs - 1) + 1; % increased number of x values
				z = imresize(z, [xs, xs]); % interpolate: xs x ys x ps
				%	x = d.xStim; x = linspace(x(1), x(end), xs); % recalculate x values
				d.(m.z) = shiftdim(z, -1); % store
				%	d.xStim = x; % store

			case 'sat' % % desaturate response colour coding
				switch m.prep.var
					case 'dir', z = d.dirPref'; % pref. dir. (deg), [-180, 180): ls x 1
					case 'max', z = d.max.val; % max. of tuning curve (mV or Hz): ls x 1
				end
				dir = d.dir; % direction (deg): 1 x ds
				dirTun = d.dirTun.val; % direction tuning (mV or Hz): ls x ds
				z = m.saturate(z, dir, dirTun, m.prep.var); % colour map: ls x 3
				d.(m.z) = shiftdim(z, -1); % RGB array: ls x 3
				d.Properties.CustomProperties.RespDim = {'', 'loc', 'colour'}; % update
			case 'scatter' % prepare for scatter plot, with one row for each stage
				d.x = struct('val', d.prop.loc(:, 1)); % horizontal loc. (deg): ls x 1
				d.y = struct('val', d.prop.loc(:, 2)); % vertical loc. (deg): ls x 1
				if isfield(d.prop, 'weight') && ~ isempty(d.prop.weight) % set size
					diam = max(1, m.prep.diam * d.prop.weight.betaEx);
						% synaptic weight, indicated by marker size
					d.size = struct('val', diam); % add it as a variable: 1 x ls
				else, d.size = struct('val', []); % default
				end
				switch d.stage % set colours
					case {'betaEx', 'gangOn', 'genOn'} % on-centre or excitatory
						d.col = struct('val', 'r'); % red for excitatory
					case {'betaIn', 'gangOff', 'genOff'} % off-centre or inhibitory
						d.col = struct('val', 'b'); % blue for inhibitory
					case 'cone' % cones
						c = eye(3); % RGB colour map for L, M, and S cones: 3 x 3
						j = d.prop.type; % cell type: ls x 1
						d.col = struct('val', c(j, :)); % marker colours: ls x 3
				end
				if isfield(m.prep, 'cone') && m.prep.cone && ...
					ismember(d.stage, {'genOff', 'genOn'}) % private line
					c = eye(3); % RGB colour map for L, M, and S cones: 3 x 3
					j = m.p.cone.type; % cone type: ls x 1
					d.col = struct('val', c(j, :)); % marker colours: ls x 3
				end
			case 'unwrap' % unwrap phase
				r = d.resp; % response: 1 x vs (vs can be multidimensional)
				dimU = m.prep.unwrap; % response dimension to unwrap
				dim = d.Properties.CustomProperties.RespDim; % response dimensions
				[~, i] = ismember(dimU, dim); % index of dimension to unwrap
				r = pi * r / 180; % convert to radians
				r = unwrap(r, [], i); % unwrap
				d.resp = 180 * r / pi; % convert back to degrees
		end
	end

function [d, m] = radius(d, m) % calculate mechanism radius

	% Calculate and store radius
	f = d.freqS; % spatial frequency (cycles/deg): 1 x fs
	r = shiftdim(d.resp, 3); % spatial frequency response (mV): 1 x fs
	fC = interp1(r, f, r(1) * exp(-1)); % interpolate for 1/e point (cycles/deg)
	rad = 1 / (pi * fC); % radius (deg)
	d.ecc = m.p.ecc; d.radius = rad; % store
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Radius (deg)';

	% Append file
	dC = d; % current data
	file = m.save.name; % name of file to which current data will be appended
	if exist(file, 'file') % there are existing data
		load(file, 'd'); % load the existing data
		d = [d; dC]; % append current data
	end

function [d, m] = resp(d, m) % calculate the model time course

% Inputs:
%		d is a single-row table containing synaptic modulation factors
%		m.resp.par, where par is a model parameter: matrix of parameter values,
%			with one row per value, e.g. m.resp.cont = [1, 1, 0; 0, 0, 1]
%		m.resp.stage: names of stages to return: char or string
%	Output:
%		a single-row table in which the response variable is:
%			a double array for a single neuronal array:
%				response (mV or Hz): ts x ls x ss x vs where vs can be multi-dim.
%			a structure array for multiple neuronal arrays (the array cannot be double
%				because neuronal arrays differ in the number of cell locations):
%				resp.resp (mV or Hz): ts x ls x vs where vs can be multi-dim.

	%	Initialise
	m = m.setWeight(d, m); % set the synaptic weights
	[name, val, vals] = m.getVal(m, 'resp'); % model parameter names and values
	vs = prod(vals); % total number of values

	%	Set the output domain and the temporal values
	if ~ isfield(m.resp, 'domain') % m.resp.domain is not set
		m.resp.domain = 'freq'; % default
	end
	switch m.resp.domain % set temporal variable
		case 'freq', t = m.p.f; % frequency (Hz): 1 x ts
		case 'time', t = m.p.t; % time (s): 1 x ts
	end
	ts = length(t); % number of frequencies or times
	
	% Initialise the response array
	i = matches(string({m.p.cell.type}), string(m.resp.stage));
		% indices of required stages in list of all stages
	s = {m.p.cell(i).type}; % stage names: 1 x ss, char
	ss = length(s); % number of stages
	a = {m.p.cell(i).array}; % array names: 1 x ss, char
	l = zeros(1, ss); % number of locations in each stage: 1 x ss
	for i = 1: ss % loop over stages
		loc = shiftdim(m.p.(a{i}).loc, -1); % neuronal locations (deg): 1 x ls x 2
		l(i) = size(loc, 2); % number of locations in stage
	end
	ls = sum(l); % total number of locations
	r = zeros(ts, ls, vs); % storage for response (mV or Hz): ts x ls x vs

	% Calculate and store the responses
	for v = 1: vs % loop over variable values

		% Set the model parameter values, and the temporal parameters
		m = m.setVal(m, name, val, v); % set values of model parameters
		m = m.samp(m); % recalculate temporal parameters

		% Calculate the response
		switch m.p.solver % choose solver
			case 'solveF' % solve for fundamental response
				p = m.solveF(m); % potential (mV), all stages in model: 1 x ss, struct
				dom = "freq"; % frequency domain response
			case 'solveT' % solve in time domain
				p = m.solveT(m); % potential (mV), all stages in model: 1 x ss, struct
				dom = "time"; % temporal domain response
		end

		%	Store the response
		i = matches({p.type}, s); % indices of stages to keep
		p = p(i); % keep required stages: 1 x ss, struct
		p = cat(1, p.resp); % concatenate potentials (mV): ls x ts, double
		p = p .'; % convenient format for temporal processing: ts x ls
		r(:, :, v) = p; % potential (mV): ts x ls x vs

	end

	% Convert to impulse rate if required
	if m.p.act % convert to impulse rate
		if dom == "freq" % potential is in the frequency domain
			r = m.ifftReal(r); % convert to temporal domain
			dom = "time"; % temporal domain response
		end
		r = m.p.kRect * max(0, r); % impulse rate (Hz)
	end

	% Set output response domain
	switch m.resp.domain % output frequency or temporal response?
		case 'freq' % frequency response
			if dom == "time", r = fft(r); end % convert temporal response
		case 'time' % temporal response
			if dom == "freq", r = ifft(r); end % convert frequency response
	end
	
	% Store variables apart from response in table
	ecc = m.p.ecc; width = m.p.wid; % eccentricity (deg), visual field width (deg)
	if isempty(d), phase = nan; cycle = nan;
	else, phase = d.phase; cycle = d.cycle; end % development phase and cycle
	d = table(ecc, width, phase, cycle, t); % store
	d.loc = loc; d.stage = s; % add location, stage
	if m.p.stimS == "image", val{1} = m.resp.index; end % replace image with index
	for i = 1: length(name) % loop over parameters
		nameC = name{i}; % current parameter
		d.(nameC) = shiftdim(val{i}, -1); % store values in table
	end

	%	Store response in table
	if isscalar(unique(a)) % single neuronal array
		d.resp = reshape(r, [1, ts, l(1), ss, vals]); % store potential:
			% 1 x ts x ls x ss x vs, where vs can be multi-dimensional
		dim = [{'', m.resp.domain, 'loc', 'stage'}, name]; % response dimensions
	else % multiple neuronal arrays: store response as structure
		iEnd = 0; locC = cell(1, ss); rC = cell(1, ss); % storage for loop
		for i = 1: ss % loop over stages
			locC{i} = m.p.(a{i}).loc; % locations (deg): ls x 2
			iStart = iEnd + 1; iEnd = iEnd + l(i); % start and end indices of location
			rC{i} = r(:, iStart: iEnd, :); % response (mV or Hz): ts x ls x vs
		end
		d.resp = struct('stage', s, 'loc', locC, 'resp', rC); % 1 x ss, struct
		dim = [{'', m.resp.domain, 'loc'}, name]; % response dimensions
	end
	
	% Set table properties
	d.Properties.VariableNames{'t'} = m.resp.domain; % make the name meaningful
	d.Properties.Description = m.p.project; % research project
	if ~ m.p.act % response quantity is generator potential
		d.Properties.VariableDescriptions{'resp'} = 'Generator potential (mV)';
	else
		d.Properties.VariableDescriptions{'resp'} = 'Impulse rate (Hz)';
	end
	d = addprop(d, 'RespDim', 'table'); % add property: response dimensions
	d.Properties.CustomProperties.RespDim = dim; % response dimensions

function [d, m] = scat(d, m) % prepare scatter plot

	d = prim.scatFunc(d, m); % prepare scatter plot

function [d, m] = spec(d, m) % calculate power spectrum of orientation map

	% Calculate orientation vector
	%	xs = size(f, 1); % number of x values
	%	xs = 2 * floor(.5 * xs); % make it even
	%	f = f(1: xs, 1: xs); % trim to make the number of x values even
	%	f = 2 * pi * f / 180; % cyclic orientation (rad): xs x ys
	%	f = exp(1i * f); % % orientation field: xs x ys
	dir = d.dir; % stimulus direction (deg), [-180, 180): 1 x ds
	r = shiftdim(d.dirTun, 1); % direction tuning (Hz): ls x ds
	ind = m.osi(dir, r)'; % orientation selectivity index: 1 x ls
	dirPref = d.dirPref; % pref. direction (deg), [-180, 180): 1 x ls
	o = m.dir2orient(dirPref); % orient. vector, length 1, variable angle: 1 x ls
	o = ind .* o; % set length of vector equal to orientation sel. index: 1 x ls

	%	Calculate sprectrum
	xs = sqrt(length(o)); % number of x locations
	o = reshape(o, [xs, xs]); % spatial array for fft2: xs x ys
	o = o .'; % make it an image: ys x xs
	p = fft2(o); % Fourier transform: xs x ys
	p = abs(p); % amplitude of Fourer transform: xs x ys
	p = p .^ 2; % power: xs x ys
	p = fftshift(p); % shift zero frequency to centre: xs x ys

	% Calculate frequency
	fs = size(p, 1); % number of frequencies
	n = .5 * fs; % number of non-negative frequencies
	i = linspace(- n, n, fs + 1); % sample indices: 1 x (fs + 1)
	i = i(1: end - 1); % open-ended interval: 1 x fs
	f = 1 / d.width; % fundamental frequency (cycles/deg)
	f = i * f; % frequencies (cycles/deg): 1 x fs

	%	Store
	d.freq = f; % frequencies (cycles/deg): 1 x fs
	d.spec = shiftdim(p, -1); % power spectrum: 1 x fs x fs

function [d, m] = stage(d, m) % list properties of specified stages

%	Assumptions:
%		input table is one row of the modulation file
%		m.stage.name = stage names: cell

	%	Construct structure of stages
	name = string(m.stage.name); % names of stages: string
	name = name(:); % make it a column vector: ss x 1
	ss = length(name); % number of stages
	p = cell(ss, 1); % empty array: ss x 1, cell
	p = struct('name', p); % stage properties: ss x 1, struct
	for i = 1: ss % loop over stages

		%	Assemble data for current stage
		nameC = name(i); % name of current stage
		j = ismember({m.p.cell.type}, nameC); % index in list of stages
		s = m.p.cell(j); % structure for current stage
		a = s.array; % name of array to which stage belongs

		%	Store properties for current stage
		p(i).name = nameC; % stage name: string
		p(i).loc = m.p.(a).loc; % cell locations (deg): ls x 2
		p(i).type = m.p.(a).type; % cell type

	end

	%	Store
	d = repmat(d, [ss, 1]); % one row for each stage
	d.stage = name; % stage name
	d.prop = p; % stage properties

function [d, m] = stim(d, m) % calculate the stimulus

	switch m.p.stimS % spatial stimulus
		case 'grating'

			%	Initialise
			t = m.p.t; % frame times: 1 x ts
			ts = length(t); % number of frames
			loc = m.listLoc(m); % make a list of grid locations: ls x 2

			%	Make the grating
			s = m.calStim(t, loc, m); % stimulus: ts x ls x cs
			x = m.p.x; % x values: 1 x xs
			xs = length(x); % number of x values
			s = reshape(s, ts, xs, xs, []); % reshape list to map: ts x xs x ys x cs
			s = permute(s, [3, 2, 4, 1]); % permute map to graph: ys x xs x cs x ts
			s = .5 * (1 + s); % change range from [-1, 1] to [0, 1] (RGB-units)

		case 'image'

			%	Initialise
			x = .5 * m.p.wid * [-1, 1]; % x limits
			t = 0; % time (s)

			%	Read image
			if isfield(m.stim, 'image'), i = m.stim.image; % index of image
			else, i = 1; end % first image only
			s = d.image(i).im; % select image: ys x xs x cs, uint8
			s = flip(s); % flip from image to graph: ys x xs x cs, uint8
			
	end

	%	Store
	s = shiftdim(s, -1); % prepare to store: 1 x ys x xs x cs x ts
	d = table(m.p.wid, x, x, t, s, 'variableNames', ...
		{'width', 'y', 'x', 'time', 'image'}); % store
	d.Properties.VariableDescriptions{'width'} = 'Visual field width (deg)';
	d = addprop(d, 'RespDim', 'table'); % add property: response dimensions
	d.Properties.CustomProperties.RespDim = ...
		{'', 'y', 'x', 'colour', 'time'}; % response dimensions

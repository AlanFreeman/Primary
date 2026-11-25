function [d, m] = index(d, m)

% Calculate response indices
%
%	Inputs:
%		some or all contrasts:
%			M-specific, equiluminant, L-specific, achromatic
%		outputs from task "field" or "max"
%	Output:
%		d.resp = structure with one field for each index, e.g. coi, soi
%	Index naming:
%		"abcD", e.g. "osiE", where:
%			a = domain/stimulus:
%				c = cone
%				e = equiluminance
%			b = analysis:
%				m = response magnitude
%				o = opponency
%				t = tuning
%			c = normalisation:
%				i = normalised index
%				u = unnormalised
%			D = optional moderator:
%				A = achromatic
%				E = equiluminant
%				R = receptive field

	%	Calculate indices
	if isfield(m.index, "name") % is an index required?
		n = string(m.index.name); % indices to calculate: 1 x is, string
		index = calIndex(d, m, n); % calculate indices
	else, return % nothing to do
	end

	%	Store
	d = d(1, :); % keep just the first row
	d.resp = struct; % initialise as structure
	for nC = n % loop over indices
		d.resp.(nC).val = index.(nC); % chromatic opponency index: ls x vs
		d.resp.(nC).dim = {'loc'}; % response dimensions *** fix ***
	end

function index = calIndex(d, m, name) % calculate response indices

%	Input:
%		name = names of indices to calculate, e.g. ["coi", "soi"]: string
%	Output:
%		index = structure with index given by field name, e.g. index.coi

	%	Assign stimulus contrasts to data file rows
	cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1]; % contrast
	stim = ["M", "equi", "L", "ach"]; % stimulus types
	for i = 1: size(d, 1) % loop over data file rows
		j = ismember(cont, d.cont(i, :), "rows"); % which stimulus type?
		d.Properties.RowNames(i) = stim(j); % assign type
	end

	%	Initialise for cone opponency index based on receptive field
	if matches("coi", name) % index required?
		iFile = d.file(1); % serial number of modulation factors file
		file = "Field " + iFile + ".mat"; % name of receptive field file
		s = string(d.stage(1)); % response stage
		if s ~= "betaEx" % not beta excitatory
			file = "Field " + s + " " + iFile + ".mat"; % name of r.f. file
		end
		s = load(file, 'd'); dR = s.d; % load receptive field data
		for i = 1: size(dR, 1) % loop over data file rows
			j = ismember(cont, dR.cont(i, :), "rows"); % which stimulus type?
			dR.Properties.RowNames(i) = stim(j); % assign type
		end
		dR = prim.crop(dR, m); % crop response to same size as that in table d
	end

	%	Initialise for structural indices
	locCone = m.p.cone.loc; % cone location (deg): cs x 2
	locCort = shiftdim(d.loc(1, :, :), 1); % cortical cell location (deg): ls x 2

	%	Prepare for composite indices
	nA = cell(1); i = 0; % additional indices, number of additional indices
	for n = name % loop over indices
		i = i + 1; % update
		switch n % composite indices
			case "bpiC" % centre-cone-specific bandpass index: Lee (12)
				nA{i} = ["bpiL", "bpiM"]; % bpiC requires these basic indices
			case "doi" % double opponency index
				nA{i} = ["coi", "soi"]; % doi requires these basic indices
			case "emuxsoi"
				nA{i} = ["emu", "soi"]; % emuxsoi requires these basic indices
			case "pca"
				nA{i} = m.pca.name; % pca requires these basic indices
		end
	end
	name = [name, nA{:}]; % include required indices
	name = unique(name); % don't duplicate indices

	%	Calculate functional indices
	for n = name % loop over indices
		switch n % choose index
			case "ami" % maximum achromatic response index
				r = d{"ach", "resp"}.max.val; % achromatic response: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				index.(n) = r; % index: ls x 1
			case "amu" % unnormalised maximum achromatic response
				index.(n) = d{"ach", "resp"}.max.val; % ls x 1
			case "ati" % achromatic tuning index
				r = d{"ach", "resp"}.freqS.tun.val; % ach. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				index.(n) = r; % ls x fs
			case {"bpiA", "bpiL", "bpiM"} % bandpass index: Lee (12)
				s = extract(n, strlength(n)); % stimulus: last character in name
				if s == "A", s = "ach"; end % achromatic stimulus
				rMax = d{s, "resp"}.max.val; % max. L-specific response: ls x 1
				r0 = d{s, "resp"}.freqS.tun.val; % achromatic tuning: ls x fs
				r0 = r0(:, 1); % response at 0 cycles/deg: ls x 1
				index.(n) = r0 ./ rMax; % bandpass index: ls x 1
			case "cdi" % centre-of-mass displacement index
				%	Displacement of cortical cell from centre of mass of L cones, M cones:
				% see prim.index
			case "coiE" % chromatic opponency index based on equiluminance response
				rE = d{"equi", "resp"}.max.val; % equiluminant response: ls x 1
				rA = d{"ach", "resp"}.max.val; % achromatic response: ls x 1
				index.(n) = rE ./ (rE + rA); % index: ls x 1
			case "coi" % chromatic opponency index based on receptive field
				rL = dR{"L", "resp"}; % L-specific response: 1 x ls x ps x 3
				rL = shiftdim(rL(:, :, :, 1), 1); % ls x ps
				rM = dR{"M", "resp"}; % M-specific response: 1 x ls x ps x 3
				rM = shiftdim(rM(:, :, :, 2), 1); % ls x ps
				index.(n) = 1 - sum(abs(rL + rM), 2) ./ sum(abs(rL) + abs(rM), 2); % ls x 1
			case "coiD" % chromatic opponency index with dominant preferred orient.

				%	Find orientation giving maximum response across both ach. and equi.
				%	*** Orientation? It looks like preferred spatial frequency is used ***
				rAT = d{"ach", "resp"}.freqS.tun.val; % achromatic s.f.r.: ls x fs
				rET = d{"equi", "resp"}.freqS.tun.val; % equiluminant s.f.r.: ls x fs
				r = [rAT, rET]; % combine ach. and equi. s.f.r.: ls x (2 * fs)
				[~, i] = max(r, [], 2); % index of maximum response: ls x 1
				[ls, fs] = size(rAT); % number of locations, frequencies
				i = mod(i, fs); % index of required resp. for ach. and equi.: ls x 1

				%	Calculate index
				rA = zeros(ls, 1); rE = rA; % allocate storage
				for j = 1: ls % loop over locations
					iC = i(j); % current index
					rA(j) = rAT(j, iC); % achromatic maximum: ls x 1
					rE(j) = rET(j, iC); % equiluminant maximum: ls x 1
				end
				index.(n) = rE ./ (rE + rA); % index: ls x 1

			case "coiS" % cone opponency index with same dir., s.f. for equi., ach.

				%	Obtain equiluminant maximum response and its indices in the full resp.
				r = d{"equi", "resp"}; % equiluminant response structure
				rE = r.max.val; % equiluminant max. response: ls x 1
				%{
				t = r.dir.pref.val; % preferred direction (deg): ls x 1
				i = knnsearch(d.dir(1, :)', t); % index of preferred direction: ls x 1
				t = r.freqS.pref.val; % preferred spatial frequency (cycles/deg): ls x 1
				j = knnsearch(d.freqS(1, :)', t); % index of pref. spat. freq.: ls x 1

				%	Obtain achromatic response, and calculate the index
				ls = length(t); rA = zeros(ls, 1); % allocate storage
				for k = 1: ls % loop over locations
					rA(k) = d{"ach", "resp"}.full.val(k, i(k), j(k)); % ach. response
				end
				%}
				rAS = respAchSame(d); % achromatic resp.: same dir., s.f. as equi. resp.
				index.(n) = rE ./ (rE + rAS); % index: ls x 1

			case "emi" % maximum equiluminant response index
				r = d{"equi", "resp"}.max.val; % equiluminant response: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				index.(n) = r; % index: ls x 1
			case "emu" % unnormalised maximum equiluminant response
				index.(n) = d{"equi", "resp"}.max.val; 
			case "eti" % equiluminant tuning index
				r = d{"equi", "resp"}.freqS.tun.val; % equi. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				index.(n) = r; % ls x fs
			case "etu" % unnormalised equiluminant tuning
				index.(n) = d{"equi", "resp"}.freqS.tun.val; % equi. tuning: ls x fs
			case "fsiE" % frequency selectivity index
				r = d{"equi", "resp"}.freqS.tun.val; % equiluminant s.f.r.: ls x fs
				%	freqS = d.freqS(1, :); % spatial frequency (cycles/deg): 1 x fs
				%	r = r(:, freqS < 20); % low freq. to avoid descending limb: ls x fs
				r = diff(r, [], 2); % approx. derivative: ls x fs, fs is less than orig.
				r = mean(r, 2); % mean of differences: ls x 1
				rMin = min(r); rMax = max(r); % minimum and maximum
				switch 2 % set bounds to [-1, 1]
					case 1 % single gradient
						g = 2 / (rMax - rMin); % gradient
						o = 1 - rMax * g; % offset
						r = o + g * r; % [-1, 1]
					case 2 % piecewise linear
						i = r < 0; r(i) = - r(i) / rMin; % set minimum to -1
						i = r >= 0; r(i) = r(i) / rMax; % set maximum to 1
				end
				r = abs(r); % absolute of mean: ls x 1
				index.(n) = r; % store: ls x 1
			case "lmi" % maximum L-specific response index
				r = d{"L", "resp"}.max.val; % response amplitude: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				index.(n) = r; % store: ls x 1
			case "lti" % L-specific tuning index
				r = d{"L", "resp"}.freqS.tun.val; % L-spec. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				index.(n) = r; % ls x fs
			case "mmi" % maximum M-specific response index
				r = d{"M", "resp"}.max.val; % response amplitude: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				index.(n) = r; % store: ls x 1
			case "mti" % M-specific tuning index
				r = d{"M", "resp"}.freqS.tun.val; % M-spec. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				index.(n) = r; % ls x fs
			case "mtu" % unnormalised M-specific tuning
				index.(n) = d{"M", "resp"}.freqS.tun.val; % M-specific s.f. tuning: ls x fs
			case {"osiA", "osiE"} % orientation selectivity index
				dir = d.dir(1, :); % stimulus direction (deg): 1 x ds
				switch n % which contrast?
					case "osiA", c = "ach"; % achromatic
					case "osiE", c = "equi"; % equiluminant
				end
				r = d{c, "resp"}.dir.tun.val; % dir. tuning (mV or Hz): ls x ds
				r = m.osi(dir, r); % orientation selectivity index, [0, 1]: ls x 1
				index.(n) = r; % store: ls x 1
			case "sdi" % S-cone distance index
				type = m.p.cone.type; % cone type: cs x 1
				locS = locCone(type == 3, :); % S-cone location (deg): ss x 2
				i = knnsearch(locS, locCort); % index of nearest S-cone: ls x 1
				dist = locCort - locS(i, :); % displacement from nearest S-cone: ls x 2
				dist = sqrt(sum(dist .^ 2, 2)); % dist. to nearest S-cone (deg): ls x 1
				index.(n) = dist / max(dist); % index: ls x 1
			case "soi" % spatial opponency index
				rLU = d{"L", "resp"}.freqS.tun.val; % spat. freq. resp. to L-spec. stim.
				rLU = rLU(:, 1); % response to uniform L-specific stimulus: ls x 1
				rMU = d{"M", "resp"}.freqS.tun.val; % spat. freq. resp. to M-spec. stim.
				rMU = rMU(:, 1); % response to uniform M-specific stimulus: ls x 1
				rL = d{"L", "resp"}.max.val; % response to optimum L-spec. stim.: ls x 1
				rM = d{"M", "resp"}.max.val; % response to optimum M-spec. stim.: ls x 1
				index.(n) = 1 - (rLU + rMU) ./ (rL + rM); % index: ls x 1
			case "x" % horizontal location
				x = locCort(:, 1); % horizontal location (deg): ls x 1
				index.(n) = x; % index: ls x 1
			case "y" % vertical location
				y = locCort(:, 2); % vertical location (deg): ls x 1
				index.(n) = y; % index: ls x 1
		end
	end

	%	Calculate composite indices
	for n = name % loop over indices
		switch n % choose index
			case "bpiC" % centre-cone-specific bandpass index: Lee (12)
				type = m.p.cone.type; % cone type
				locC = shiftdim(d.loc(1, :, :), 1); % cropped locations (deg): ls x 2
				ls = size(locC, 1); % number of cropped cells
				bpi = [index.bpiL, index.bpiM, nan(ls, 1)];
					% cone-specific bandpass indices: ls x 3
				ind = zeros(ls, 1); % allocate storage
				for i = 1: ls % loop over cells
					j = knnsearch(locCone, locC(i, :)); % index of current cone
					typeC = type(j); % type of current cone
					ind(i) = bpi(i, typeC); % bandpass index for stimulation of that cone
				end
				index.bpiC = ind; % store: ls x 1
			case "doi" % double opponency index
				index.(n) = index.coi .* index.soi; % product
			case "emuxsoi" % product of equiluminance magnitude and soi
				index.(n) = index.emu .* index.soi; % requires emu and soi
			case "pca" % indices for principal components analysis
				nC = string(m.pca.name); % indices providing variables for PCA
				ns = length(nC); % number of indices
				r = cell(1, ns); % allocate storage
				for i = 1: ns % loop over indices
					r{i} = index.(nC(i)); % ith index
				end
				index.(n) =[r{:}]; % concatenate
		end
	end

function rA = respAchSame(d) % achromatic resp.: same dir., s.f. as equi. resp.

	%	Obtain equiluminant maximum response and its indices in the full resp.
	r = d{"equi", "resp"}; % equiluminant response structure
	p = r.dir.pref.val; % preferred direction (deg): ls x 1
	i = knnsearch(d.dir(1, :)', p); % index of preferred direction: ls x 1
	p = r.freqS.pref.val; % preferred spatial frequency (cycles/deg): ls x 1
	j = knnsearch(d.freqS(1, :)', p); % index of pref. spat. freq.: ls x 1

	%	Obtain achromatic response, and calculate the index
	ls = length(p); rA = zeros(ls, 1); % allocate storage
	for k = 1: ls % loop over locations
		rA(k) = d{"ach", "resp"}.full.val(k, i(k), j(k)); % ach. response: ls x 1
	end
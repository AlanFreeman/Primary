function index = indexFunc(d, m, name) % calculate response indices

%	Input:
%		name = names of indices to calculate, e.g. ["coi", "soi"]: string
%	Output:
%		index = structure with index given by field name, e.g. index.coi

	%	Assign stimulus types to data file rows
	cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1]; % contrast
	stim = ["M", "equi", "L", "ach"]; % stimulus types
	for i = 1: size(d, 1) % loop over data file rows
		j = ismember(cont, d.cont(i, :), "rows"); % which stimulus type?
		d.Properties.RowNames(i) = stim(j); % assign type
	end

	%	Calculate functional indices
	for n = name % loop over indices
		switch n % choose index
			case "ami" % maximum achromatic response index
				r = d{"ach", "resp"}.max.val; % achromatic response: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				ind = r; % index: ls x 1
			case "ati" % achomatic tuning index
				r = d{"ach", "resp"}.freqS.tun.val; % ach. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				ind = r; % ls x fs
			case "coi" % chromatic opponency index
				rE = d{"equi", "resp"}.max.val; % equiluminant response: ls x 1
				rA = d{"ach", "resp"}.max.val; % achromatic response: ls x 1
				ind = rE ./ (rE + rA); % index: ls x 1
			case "coiDom" % chromatic opponency index with dominant preferred orient.

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
				ind = rE ./ (rE + rA); % index: ls x 1

			case "coiE" % cone opponency index with same dir., s.f. for equi., ach.

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
				ind = rE ./ (rE + rAS); % index: ls x 1

			case "emi" % maximum equiluminant response index
				r = d{"equi", "resp"}.max.val; % equiluminant response: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				ind = r; % index: ls x 1
			case "eti" % equiluminant tuning index
				r = d{"equi", "resp"}.freqS.tun.val; % equi. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				ind = r; % ls x fs
			case "fsi" % frequency selectivity index
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
				ind = r; % store: ls x 1
			case "lmi" % maximum L-specific response index
				r = d{"L", "resp"}.max.val; % response amplitude: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				ind = r; % store: ls x 1
			case "lti" % L-specific tuning index
				r = d{"L", "resp"}.freqS.tun.val; % L-spec. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				ind = r; % ls x fs
			case {"amu", "emu"} % unnomalised maximum response
				switch n % which contrast?
					case "amu", c = "ach"; % achromatic
					case "emu", c = "equi"; % equiluminant
				end
				ind = d{c, "resp"}.max.val; % maximum response: ls x 1
			case "mmi" % maximum M-specific response index
				r = d{"M", "resp"}.max.val; % response amplitude: ls x 1
				r = r / max(r); % normalised response amplitude: ls x 1
				ind = r; % store: ls x 1
			case "mti" % M-specific tuning index
				r = d{"M", "resp"}.freqS.tun.val; % M-spec. spat. freq. tuning: ls x fs
				r = r ./ max(r, [], 2); % normalise: ls x fs
				ind = r; % ls x fs
			case {"osiA", "osiE"} % orientation selectivity index
				dir = d.dir(1, :); % stimulus direction (deg): 1 x ds
				switch n % which contrast?
					case "osiA", c = "ach"; % achromatic
					case "osiE", c = "equi"; % equiluminant
				end
				r = d{c, "resp"}.dir.tun.val; % dir. tuning (mV or Hz): ls x ds
				r = m.osi(dir, r); % orientation selectivity index, [0, 1]: ls x 1
				ind = r; % store: ls x 1
			case "soi" % spatial opponency index
				rLU = d{"L", "resp"}.freqS.tun.val; % spat. freq. resp. to L-spec. stim.
				rLU = rLU(:, 1); % response to uniform L-specific stimulus: ls x 1
				rMU = d{"M", "resp"}.freqS.tun.val; % spat. freq. resp. to M-spec. stim.
				rMU = rMU(:, 1); % response to uniform M-specific stimulus: ls x 1
				rL = d{"L", "resp"}.max.val; % response to optimum L-spec. stim.: ls x 1
				rM = d{"M", "resp"}.max.val; % response to optimum M-spec. stim.: ls x 1
				ind = 1 - (rLU + rMU) ./ (rL + rM); % index: ls x 1
			case "x" % horizontal location
				x = shiftdim(d.loc(1, :, 1), 1); % horizontal location (deg): ls x 1
				ind = x; % index: ls x 1
			case "y" % vertical location
				y = shiftdim(d.loc(1, :, 2), 1); % vertical location (deg): ls x 1
				ind = y; % index: ls x 1
		end
		index.(n) = ind; % ls x ps
	end

	%	Calculate composite indices
	for n = name % loop over indices
		switch n % choose index
			case "doi" % double opponency index
				index.(n) = index.coi .* index.soi; % requires coi and soi
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

	%	Calculate structural indices
	for n = name % loop over indices
		switch n % choose index
			case {"cdi", "sdi"}
				index.(n) = prim.indexStruct(d, m, n); % ls x ps
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

function m = setLit(m) % set metadata calculated from the literature

	% Files:
	%		'Functional eccentricity.mat'

	% Parameters calculated from the literature, e = functional eccentricity (deg)
	%	m.p.densBeta = [178380, 32720]; % dens. of 4CBeta exc., inh. cells (mm^-3)
		%	O'Kusky (82), Fitzpatrick (87)
	m.p.actRest = [23.6, 23.6]; % resting gang. cell impulse rate (Hz): [off, on]
		% Troy (94)
	m.p.densBeta = 3e4; % density of 4CBeta neurons (cells/mm^2): O'Kusky (82)
	%	m.p.densGangCoef = [1.3497, -7.8485, 7.7392, 5.6107]; % g.c. density
		%	Wässle (89)
		% ganglion cell density (deg^-2) = exp(polyval(c, log10(eAnat)))
	m.p.kGangDev = .1241; % s.d. of gang. cell nearest-neighbour dist., Dacey (93)
	m.p.kGangGen = .47; % attenuation from g.c. to l.g.n., Kaplan (87)
	m.p.kRect = 7.2; % rectification constant (Hz/mV), Carandini (00)
	m.p.kSens = 17.8; % gang. cell contrast sens'y (mV/cont.-unit), Croner (95)
	m.p.kSensBeta = 23; % exc. cell cont. sens. (Hz/contrast-unit), Edwards (95);
		% an estimate from layer 2/3: find a better estimate
	m.p.kSur = 1.54; % surround gain, Croner (95)
	m.p.magRet = 4.73; % retinal magnification factor (deg/mm), Perry (85)
	%	m.p.potRest = 5.7; % g.c. resting pot. (mV), Kaplan (87): 41.1 / m.p.kRect
	%	m.p.radFix = .1; % radius over which to sum image background (deg):
		% an estimate from fixation eye movements: Martinez-Conde (13)
	m.p.radFix = .127; % radius over which to sum image background (deg):
		% an estimate from fixation eye movements: Skavenski (75)
	m.p.radBetaMm = .1421; % 4CBeta convergence radius (mm)
	%	m.p.radGangCoef = [-1.3583e-06, .00014598, .0012385, .012601]; % g.c. radius
		%	Watanabe (89)
		% ganglion cell dendritic radius (deg) = polyval(c, e)
	m.p.radHorCoef = [-4.5667e-05, .0030709, .029977, -3.261]; % h.c. r.f. radius
		% Wässle (89), Packer (02)
		% horizontal cell receptive field radius (deg) = exp(polyval(c, e))
	m.p.radInExMm = .154; % inhibitory to excitatory convergence radius (mm):
		% mean mouse somatosensory cortex, Packer (11)
	%	m.p.radOptCoef = [1.4241e-05, 0, .011663]; % point spread function radius
	m.p.radOptCoef = [1.8988e-05, 0, .015551]; % point spread function radius
		% Navarro (93), Godat (22)
		% point spread function radius (deg) = polyval(c, e)
	m.p.ratBetaIn = .143; % ratio of inhib. to all beta neurons, Fitzpatrick (87)
	m.p.ratCone = [0.4439, 0.4310, 0.1251]; % cone ratio: [L, M, S] / (L + M + S]
		%	Munds (22)
	m.p.ratGang = .835; % ratio of midget to all ganglion cells %	Peng (19)
	%	m.p.ratOpt = 1.5; % ratio of PSF radius, macaque/human, Harwerth (85)
	m.p.ratOpt = 2; % ratio of PSF radius, macaque/human, match Williams (81)
	%	m.p.ratSign = .581; % ratio of off- to all midget ganglion cells, Peng (19)
	m.p.ratSign = .52; % ratio of off- to all midget ganglion cells, pragmatic
	m.p.tau = .0123; % time constant (s), Ringach (97)
	m.p.xyz2lms = [ % conversion from XYZ colour space to LMS: Stockman (23)
		0.2106, 0.8551, -0.0397
		-0.4171, 1.1773, 0.0786
		0,      0,      0.5168];

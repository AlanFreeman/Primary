function m = setDef(m) % set default model parameters

	% Stimulus parameters
	m.p.contMag = .3; % contrast magnitude
	m.p.cont = m.p.contMag * [1, 1, 1]; % cone contrast, [cL, cM, cS]
	m.p.dir = 0; % stimulus direction (deg anticlockwise from rightward)
	m.p.dur = .05; % stimulus duration (s)
	m.p.freqS = 4; % spatial frequency (cycles/deg)
	m.p.freqT = 2; % temporal frequency (cycles/s)
	m.p.locS = 0; % displacement in stimulus direction (deg)
	m.p.phaseS = 0; % spatial phase (deg)
	m.p.stimS = 'grating'; % spatial stimulus: grating, image
	m.p.stimT = 'drift'; % temporal stimulus: drift, pulse

	%	Sampling parameters
	m.p.time = 1 / m.p.freqT; % simulation time (s) for cyclic stimulus
	m.p.ts = 128; % number of sample times and temporal frequencies: make it even
	m.p.xs = 51; % number of sample locations
	
	% Structural parameters
	m.p.array = {'cone', 'gangOff', 'gangOn', 'betaEx', 'betaIn'}; % neuron arrays
	m.p.betaEx.stage = {'betaEx', 'sumGen', 'sumIn'}; % 4CBeta excitatory stages
	m.p.betaIn.stage = {'betaIn'}; % 4CBeta inhibitory stage
	%	m.p.cone.stage = {'cone', 'hor', 'back', 'biOff', 'biOn'}; % cone stages
	m.p.cone.stage = {'coneInt', 'cone', 'hor', 'back', 'biOff', 'biOn'}; % cone
	m.p.gangOff.stage = {'gangOff', 'genOff'}; % off-centre stages
	m.p.gangOn.stage = {'gangOn', 'genOn'}; % off-centre stages
	m.p.synHebb = ["gen", "betaEx"; "gen", "betaIn"]; % Hebbian synapses
	m.p.synIn = ["betaIn", "betaEx"]; % syn. with spike-timing-independent plast.

	%	Visual field parameters
	m.p.ecc = 1; % functional eccentricity (deg)
	m.p.wid = .2; % visual field width and height (deg)
	
	%	Development parameters
	m.p.cycle = 1: 3; % development cycles written to data file
	m.p.kIncIn = .0015; % increment in inhibitory-excitatory gain per cycle
	m.p.kIncMod = 1.5; % maximum change of modulation factor per cycle

	% Pragmatic parameters
	m.p.coneStages = 2; % number of cone stages
	m.p.kAdapt = .975; % gain multiplier of L-cones when adapted
	m.p.kDensBeta = .01785; % multiple of 4CBeta cell density: reduce comp. time
	%	m.p.kDevGain = .1; % standard dev. of gain from geniculate to layer 4CBeta
	m.p.kGainOff = 1.025; % off-centre channel gain multiplier
	m.p.kGenBeta = 2; % gain from geniculate to layer 4CBeta
	%	m.p.kInEx = 1; % gain from inhibitory to excitatory cells
	%	m.p.ratSign = .52; % ratio of off- to all midgets; empirical value is .581

	% Computational flags and parameters
	m.p.act = 1; % signal type: 0 for potential, 1 for impulse rate
	m.p.adapt = 1; % L-cone adaptation: 0 for none, 1 for some
	m.p.align = 0; % alignment of gang. cells with cones: 0 for random, 1 align
	m.p.back = 1; % feedback loop: 0 for open, 1 for closed
	m.p.file = ''; % data file
	m.p.parallel = 0; % parallel pool: 0 for no, 1 for yes
	m.p.project = 'Prim'; % default data folder: userpath/MATLAB/Data/m.p.project
	m.p.seed = 0; % seed for random number generator
	m.p.solver = 'solveF'; % solver function: solveF or solveT
	m.p.varyGain = 1; % vary off-centre gain: 0 for no, 1 for yes
	m.p.warning = 1; % warning for large imag. components? 0 for no, 1 for yes

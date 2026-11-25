function p = solveT(m) % solve the model equations in the time domain

	% Set initial values for potentials
	for c = m.p.stage % loop over stages
		s = c.name; % current stage
		a = c.array; % current location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(s) = zeros(ls, 1); % default resting potential (mV): ls x 1
	end
	pRest = m.p.potRest; % ganglion cell resting potential (mV)
	p.gangOff(:) = pRest(1); p.gangOn(:) = pRest(2); % ganglion cells
	pRest = m.p.kGangGen * pRest; % geniculate resting potential
	p.genOff(:) = pRest(1); p.genOn(:) = pRest(2); % geniculate potential
	p.gen = [p.genOff; p.genOn]; % combined off- and on-centre
	pRest = m.p.kGenBeta * p.gen; % geniculate output
	pEx = m.p.w_gen_betaEx * pRest; % excitatory input
	p.betaIn(:) = m.p.w_gen_betaIn * pRest; % inhibitory potential
	p.betaEx = pEx - m.p.kInEx * m.p.w_betaIn_betaEx * p.betaIn; % excitatory pot.
	pI = struct2cell(p); % make it a cell array
	pI = vertcat(pI{:}); % concatenate across stages: ls x 1, all locations

	% Numerically integrate the equations
	fun = @(t, p) prim.deriv(t, p, m); % function to calculate derivatives
	t = m.p.t; % sampling times (s): 1 x ts
	[~, pot] = ode45(fun, t, pI); % potential (mV): ts x ls
	
	% Unpack potential
	i = 0; % current column in potential
	for s = m.p.stage % loop over stages
		sC = s.name; % current stage
		ls = length(p.(sC)); % number of locations
		pC = pot(:, i + (1: ls)); % store potentials in named array: ts x ls
		p.(sC) = pC'; % transpose to standard dimension order: ls x ts
		i = ls + i; % update column
	end

	% Add unintegrated signals
	p.back = m.p.kSur * m.p.wHorCone * p.hor; % feedback (mV):
		%	(ls x ls) x (ls x ts) = ls x ts
	p.gen = [p.genOff; p.genOn]; % combine off- and on-inputs
	pExc = p.gen; % excitatory input
	pExc = max(0, pExc); % presynaptic potential, rectified by the geniculate
	p.sumGen = m.p.kGenBeta * m.p.w_gen_betaEx * pExc; % summed excit. input
	pInh = max(0, p.betaIn); % inhibitory potential, rectified by inh. cells
	p.sumIn = m.p.kInEx * m.p.w_betaIn_betaEx * pInh; % summed inhibitory input

	% Rearrange potential to match old format
	s = m.p.stage; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).name; % current stage
		s(i).resp = p.(sC); % add response: ls x ts
	end
	p = s; % structure, including responses: 1 x zs
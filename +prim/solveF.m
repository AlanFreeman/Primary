function p = solveF(m) % solve the model equations in the frequency domain

	% Initialise
	f = 2 * pi * m.p.f; % fft frequencies (radians/s): 1 x fs
	T = 1 ./ (1 + 1i * m.p.tau * f); % temporal attenuation (radians): 1 x fs
	D = - m.p.kSens * m.driveF(m); % drive (mV): ls x fs
	w = m.p.back * m.p.kSur * m.p.wHorCone * m.p.wConeHor; % spat. filter: ls x ls
	ls = size(w, 1); % number of cones
	fs = length(f); % number of frequencies

	%	Initialise for for-loop
	if any(D(:, 3: end - 1), "all") % speed up calculation for drifting grating
		range = 1: fs; % all frequencies
	else
		range = [1, 2, fs]; % mean and fundamental components
	end
	pC = zeros(ls, fs); % allocate storage: ls x fs

	% Calculate cone signal
	switch m.p.pedicle % feedback location
		case 0 % feedback at outer sement
			p.coneInt = pC; % intermediate cone signal; no need to calculate
			n = 3; % power of temporal filter in surround
		case 1 % feedback at pedicle, not outer segment
			p.coneInt = T .* D; % intermediate cone signal
			n = 2; % power of temporal filter in surround
	end
	for j = range % loop over frequencies
		TC = T(j); % current temporal filter (radians): 1 x 1
		DC = D(:, j); % current drive (mV): ls x 1
		pC(:, j) = (eye(ls) + TC .^ n * w) \ (TC .^ 2 .* DC);
			% cone signal (mV): ls x 1
	end
	p.cone = pC; % cone signal (mV): ls x fs
	
	% Calculate remaining cone array signals
	p.hor = T .* (m.p.wConeHor * p.cone); % horizontal cell signal (mV): ls x fs
	p.back = m.p.kSur * m.p.wHorCone * p.hor; % feedback signal (mV): ls x fs
	p.biOff = T .* p.cone; % off-bipolar signal (mV): ls x fs
	p.biOn = - T .* p.cone; % on-bipolar signal (mV): ls x fs
	
	% Calculate ganglion cell array signals
	j = 0; % index for resting potential
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i; gen = "gen" + i;
			% construct weight name
		p.(gang) = T .* (m.p.(w) * p.(bi)); % ganglion cell fund. (mV): ls x fs
		j = j + 1; % update index
		p.(gang)(:, 1) = p.(gang)(:, 1) + fs * m.p.potRest(j); % rest. pot.: ls x fs
		p.(gen) = m.p.kGangGen * T .* m.rect(p.(gang), 2); % gen. pot. (mV): ls x fs
	end

	% Calculate layer 4CBeta array signals
	pGen = [p.genOff; p.genOn]; % combine off- and on-centre inputs (mV): ls x fs
	p.gen = pGen; % geniculate potential (mV): ls x fs
	pGen = m.rect(pGen, 2); % rectify input (mV): ls x fs
	pEx = m.p.kGenBeta * m.p.w_gen_betaEx * pGen; % summed exc. input (mV): ls x fs
	p.sumGen = pEx; % store summed excitation (mV): ls x fs
	pIn = m.p.kGenBeta * m.p.w_gen_betaIn * pGen;
	p.betaIn = T .* pIn; % inhibitory cell response (mV): ls x fs
	pIn = m.p.kInEx * m.p.w_betaIn_betaEx * m.rect(p.betaIn, 2); % summed: ls x fs
	p.sumIn = pIn; % store summed inhibition (mV): ls x fs
	p.betaEx = T .* (pEx - pIn); % excitatory cell response (mV): ls x fs

	% Add response to each stage
	s = m.p.stage; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).name; % current stage
		s(i).resp = p.(sC); % add response
	end
	p = s; % 1 x zs
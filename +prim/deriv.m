function dpot = deriv(t, pot, m)

% Calculate derivatives of model equations

	% Sort potential into cell types
	i = 0; % current row in potential vector
	for c = m.p.stage % loop over stages
		type = c.name; % stage name
		a = c.array; % location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(type) = pot(i + (1: ls)); % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end

	% Calculate derivatives for cone array
	rate = 1 / m.p.tau; % reciprocal of time constant (s^-1)
	drive = - m.p.kSens * m.driveT(t, m)'; % cone generator potential (mV): ls x 1
	back = m.p.back * m.p.kSur * m.p.wHorCone * p.hor; % f'back pot. (mV): ls x 1
	switch m.p.pedicle % feedback location
		case 0 % outer segment
			dp.coneInt = rate * (drive - back - p.coneInt);
				% intermediate cone potential (mV/s): ls x 1
			dp.cone = rate * (p.coneInt - p.cone); % cone potential (mV/s): ls x 1
		case 1 % pedicle
			dp.coneInt = rate * (drive - p.coneInt);
				% intermediate cone potential (mV/s): ls x 1
			dp.cone = rate * (p.coneInt - back - p.cone); % cone pot. (mV/s): ls x 1
	end
	dp.hor = rate * (m.p.wConeHor * p.cone - p.hor); % hor. cell (mV/s): ls x 1
	dp.back = 0 * dp.hor; % feedback (mV/s): ls x 1
	dp.biOff = rate * (p.cone - p.biOff); % off-bipolar (mV/s): ls x 1
	dp.biOn = rate * (- p.cone - p.biOn); % on-bipolar (mV/s): ls x 1
	
	% Calculate derivatives for ganglion cell arrays
	j = 0; % index for resting potential
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i; gen = "gen" + i;
			% construct names
		j = j + 1; % update index
		dp.(gang) = rate * (m.p.(w) * p.(bi) + m.p.potRest(j) - p.(gang));
			% ganglion cell potential (mV): ls x 1
		pC = max(0, p.(gang)); % presynaptic pot. at geniculate (mV)
		dp.(gen) = rate * (m.p.kGangGen * pC - p.(gen)); % gen. pot. (mV): ls x 1
	end

	% Calculate derivatives for layer 4CBeta
	pGen = [p.genOff; p.genOn]; % combine off- and on-inputs: ls x 1
	dp.gen = 0 * pGen; % place-holder: ls x 1
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
	for c = m.p.stage % loop over stages
		type = c.name; % name of current stage
		dpC = dp.(type); ls = size(dpC, 1); % derivatives, number of derivatives
		dpot(i + (1: ls)) = dpC; % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end
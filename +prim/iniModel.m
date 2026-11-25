function m = iniModel(m) % initialise model

	% Initialise randomisation and warnings
	rng(m.p.seed); % make the results reproducible
	i = 'runSup:LargeImaginary'; % warning about large imag. components in ifft
	if m.p.warning, warning('on', i); % define warning, turn it on
	else, warning('off', i); end % define warning, turn it off
	
	%	Initialise stages, synapses, and samples
	m = prim.iniStruct(m); % calculate radii, arrays, and convergence functions
	m = m.setSyn(m); % create list of all plastic synapse
	m = m.samp(m); % create spatial and temporal samples

	%	Initialise computed parameters from the literature
	m.p.potRest = m.p.actRest / m.p.kRect;
		% ganglion cell resting potential (mV), [off, on]: 2 x 1
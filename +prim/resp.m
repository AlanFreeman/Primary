function [d, m] = resp(d, m) % calculate the model time course

% Inputs:
%		d is a single-row table containing synaptic modulation factors
%		m.resp.par, where par is a model parameter: matrix of parameter values,
%			with one row per value, e.g. m.resp.cont = [1, 1, 0; 0, 0, 1]
%		m.resp.stage: names of stages to return: character vector or string
%	Output:
%		a single-row table in which the response variable is:
%			a double array for a single neuronal array:
%				response (mV or Hz): ts x ls x ss x vs where vs can be multi-dim.
%			a structure array for multiple neuronal arrays (the array cannot be double
%				because neuronal arrays differ in the number of cell locations):
%				resp.resp (mV or Hz): ts x ls x vs where vs can be multi-dim.

	%	Initialise
	m = m.setWeight(d, m); % set the synaptic weights
	[name, val, vals] = prim.getVal(m, "resp"); % stim. parameter names and values
	vs = prod(vals); % total number of values

	%	Set the output domain and the temporal values
	if ~ isfield(m.resp, "domain") % m.resp.domain is not set
		m.resp.domain = "freq"; % default
	end
	switch m.resp.domain % set temporal variable
		case "freq", t = m.p.f; % frequency (Hz): 1 x ts
		case "time", t = m.p.t; % time (s): 1 x ts
	end
	ts = length(t); % number of frequencies or times
	
	% Initialise the response array
	i = matches([m.p.stage.name], string(m.resp.stage));
		% indices of required stages in list of all stages
	s = [m.p.stage(i).name]; % stage names: 1 x ss, string
	ss = length(s); % number of stages
	a = [m.p.stage(i).array]; % array names: 1 x ss, string
	l = zeros(1, ss); % number of locations in each stage: 1 x ss
	for i = 1: ss % loop over stages
		loc = shiftdim(m.p.(a(i)).loc, -1); % neuronal locations (deg): 1 x ls x 2
		l(i) = size(loc, 2); % number of locations in stage
	end
	ls = sum(l); % total number of locations
	r = zeros(ts, ls, vs); % storage for response (mV or Hz): ts x ls x vs

	% Calculate and store the responses
	for v = 1: vs % loop over stimulus values

		% Set the stimulus parameter values, and the temporal parameters
		m = m.setVal(m, name, val, v); % set values of model parameters
		m = m.samp(m); % recalculate temporal parameters

		% Calculate the response
		switch m.p.solver % choose solver
			case "solveF" % solve for fundamental response
				p = prim.solveF(m); % potential (mV), all stages in model: 1 x ss, struct
				dom = "freq"; % frequency domain response
			case "solveT" % solve in time domain
				p = prim.solveT(m); % potential (mV), all stages in model: 1 x ss, struct
				dom = "time"; % temporal domain response
		end

		%	Store the response
		i = matches([p.name], s); % indices of stages to keep
		p = p(i); % keep required stages: 1 x ss, struct
		p = cat(1, p.resp); % concatenate potentials (mV): ls x ts, double
		p = p .'; % convenient format for temporal processing: ts x ls
		r(:, :, v) = p; % potential (mV): ts x ls x vs

	end

	% Convert to impulse rate if required
	if m.p.act % convert to impulse rate
		if dom == "freq" % potential is in the frequency domain
			r = prim.ifftReal(r); % convert to temporal domain
			dom = "time"; % temporal domain response
		end
		r = m.p.kRect * max(0, r); % impulse rate (Hz)
	end

	% Set output response domain
	switch m.resp.domain % output frequency or temporal response?
		case "freq" % frequency response
			if dom == "time", r = fft(r); end % convert temporal response
		case "time" % temporal response
			if dom == "freq", r = ifft(r); end % convert frequency response
	end
	
	% Store variables apart from response in table
	file = d.file; % modulation factors file serial number
	ecc = m.p.ecc; width = m.p.wid; % eccentricity (deg), visual field width (deg)
	if isempty(d), phase = nan; cycle = nan;
	else, phase = d.phase; cycle = d.cycle; end % development phase and cycle
	d = table(file, ecc, width, phase, cycle, t); % store
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
		dim = ["", m.resp.domain, "loc", "stage", name]; % response dimensions
	else % multiple neuronal arrays: store response as structure
		iEnd = 0; locC = cell(1, ss); rC = cell(1, ss); % storage for loop
		for i = 1: ss % loop over stages
			locC{i} = m.p.(a{i}).loc; % locations (deg): ls x 2
			iStart = iEnd + 1; iEnd = iEnd + l(i); % start and end indices of location
			rC{i} = r(:, iStart: iEnd, :); % response (mV or Hz): ts x ls x vs
		end
		d.resp = struct("stage", s, "loc", locC, "resp", rC); % 1 x ss, struct
		dim = ["", m.resp.domain, "loc", name]; % response dimensions
	end
	
	% Set table properties
	%	d.Properties.VariableNames("t") = m.resp.domain; % make the name meaningful
	d = renamevars(d, "t", m.resp.domain); % make the name meaningful
	d.Properties.Description = m.p.project; % research project
	if ~ m.p.act % response quantity is generator potential
		d.Properties.VariableDescriptions("resp") = "Generator potential (mV)";
	else
		d.Properties.VariableDescriptions("resp") = "Impulse rate (Hz)";
	end
	d = addprop(d, "RespDim", "table"); % add property: response dimensions
	d.Properties.CustomProperties.RespDim = dim; % response dimensions
function m = weight(m)

% Calculate synaptic convergence function attenuations

	% Calculate cone to horizontal attenuation
	locI = m.p.cone.loc; % cone locations: is x 2
	locO = locI; % horizontal cells are colocated with cones: os x 2
	r = m.p.radHor; % horizontal cell field radius
	a = m.calConv(locI, locO, r); % distance-based attenuation: os x is
	
	% Set the weight to zero for S-cone inputs, and store
	iCone = m.p.cone.type == 3; % cone is S-type: is x 1
	w = a; w(:, iCone) = 0; % H1 horizontal cells receive no S-cone input: os x is
	w = w ./ sum(w, 2); % sum of input weights is 1: os x is
	m.p.wConeHor = w; % cone to horizontal cell synaptic weights: os x is
	w = a'; w(:, iCone) = 0; % feedback connection is same as feedforward: os x is
	w = w ./ sum(w, 2); % sum of input weights is 1: os x is
	m.p.wHorCone = w; % horizontal cell to cone synaptic weights: os x is

	% Calculate bipolar to off-ganglion cell attentuation
	locI = m.p.cone.loc; % bipolar cell locations: is x 2
	locO = m.p.gangOff.loc; % ganglion cell locations: os x 2
	radGang = m.p.radGang; % ganglion cell dendritic radius
	w = m.calConv(locI, locO, radGang); % attenuation for off-cells: os x is
	w = w ./ sum(w, 2); % sum of input weights is 1: os x is
	m.p.wBiGangOff = w; % bipolar cell to off-gang. synaptic weights: os x is
	
	% Calculate bipolar to on-ganglion cell attentuation	
	locO = m.p.gangOn.loc; % ganglion cell locations: os x 2
	w = m.calConv(locI, locO, radGang); % attenuation for on-cells: os x is
	w(:, iCone) = 0; % there are no on-S midget ganglion cells: os x is
	w = w ./ sum(w, 2); % sum of input weights is 1: os x is
	if radGang == 0 % private line
		i = isnan(w); % nan results from dividing 0 by 0
		w(i) = 0; % weights are 0
	end
	m.p.wBiGangOn = w; % bipolar cell to on-gang. synaptic weights: os x is

	% Calculate geniculate to 4CBeta excitatory cell attentuation
	locI = [m.p.gangOff.loc; m.p.gangOn.loc]; % geniculate cell locations: is x 2
	locO = m.p.betaEx.loc; % 4CBeta excitatory locations: os x 2
	r = m.p.radBeta; % geniculocortical convergence radius (deg)
	a = m.calConv(locI, locO, r); % attenuation: os x is
	if radGang == 0 % private line
		a(:, iCone) = 0; % off-S geniculate cells do not project to 4CBeta
	end
	a = a ./ sum(a, 2); % sum of input attenuations is 1: os x is
	m.p.a_gen_betaEx = a; % geniculocortical attenuation: os x is

	% Calculate geniculate to 4CBeta inhibitory cell attentuation
	locO = m.p.betaIn.loc; % 4CBeta inhibitory locations: os x 2
	a = m.calConv(locI, locO, r); % attenuation: os x is
	if radGang == 0 % private line
		a(:, iCone) = 0; % off-S geniculate cells do not project to 4CBeta
	end
	a = a ./ sum(a, 2); % sum of input attenuations is 1: os x is
	m.p.a_gen_betaIn = a; % geniculocortical attenuation: os x is

	% Calculate 4CBeta inhibitory to excitatory attenuation
	locI = m.p.betaIn.loc; % 4CBeta inhibitory cell locations: is x 2
	locO = m.p.betaEx.loc; % 4CBeta excitatory cell locations: os x 2
	r = m.p.radInEx; % inhibitory-excitatory convergence radius (deg)
	a = m.calConv(locI, locO, r); % attenutation: os x is
	a = a ./ sum(a, 2); % sum of input weights is 1: os x is
	m.p.w_betaIn_betaEx = a; % inhibitory to excitatory attenuation: os x is
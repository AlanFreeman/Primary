function [d, m] = densGang(d, m) % calculate ganglion cell density v. func. ana.

	% Integrate density to find cell count
	eA = d.ecAna; % anatomical eccentricity (deg): es x 1
	dens = d.dens; % density (deg^-2): es x 1
	dn = 2 * pi * eA .* dens; % derivative of cell count (deg^-1): es x 1
	n = cumtrapz(eA, dn); % integrate to find cell count: es x 1
	
	% Convert anatomical to functional eccentricity
	eF = model.readTab("ecFunc", eA); % functional ecccentricity (deg)

	% Differentiate cell count to find density versus functional eccentricity
	[eF, i] = unique(eF); % remove repeated zero eccentricities
	es = length(eF); % number of eccentricities
	n = n(i); % remove corresponding counts
	eR = linspace(eF(1), eF(end), es); % resample eF at regular spacing: 1 x es
	n = interp1(eF, n, eR); % interpolate to obtain density: 1 x es
	dens = gradient(n, eR(2) - eR(1)) ./ (2 * pi * eR); % density is derivative

	%	Select midget cells. Obtain cone density, for low eccentricity calculations
	e = [0: .25: .75, 1: .5: 3.5, 4: 1: 9, 10: 5: 50];
		% resample eccentricity for efficient storage (deg): 1 x es
	dens = interp1(eR, dens, e); % resample density
	dens = m.p.ratGang * dens; % midget ganglion cell density (deg^-2): 1 x es
	densCone = model.readTab("densCone", e); % cone density (deg^-2): 1 x es

	%	Calculate off- and on-cell densities
	densOff = m.p.ratSign * dens; % density at ec. beyond transition interval
	densOff = model.trans(e, densCone, densOff, m); % fix densities at low ec.
	densOn = (1 - m.p.ratSign) * dens; % density at ec. beyond transition interval
	%	r = sum(m.p.ratCone(1: 2)); % fraction of cones that are L or M
	%	densOn = model.trans(e, r * densCone, densOn, m); % fix densities at low ec.
	densOn = model.trans(e, densCone, densOn, m); % fix densities at low ec.
	
	% Store
	ec = e'; densOff = densOff'; densOn = densOn'; % transpose: es x 1
	d = table(ec, densOff, densOn); % store in table: es x 3
	d.Properties.VariableDescriptions("ec") = "Functional eccentricity (deg)";
	d.Properties.VariableDescriptions("densOff") = ...
		"Off-centre midget ganglion cell density (deg^-^2)";
	d.Properties.VariableDescriptions("densOn") = ...
		"On-centre midget ganglion cell density (deg^-^2)";
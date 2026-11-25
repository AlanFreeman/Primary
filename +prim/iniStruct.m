function m = iniStruct(m) % initialise model structure

	%	Create a structure of stages
	list = stageList; % list of stages and arrays: ss x 1, string
	ss = size(list, 1); % number of stages
	s = struct; % specify class
	for i = 1: ss % loop over stages
		s(i).name = list(i, 1); % stage name
		s(i).array = list(i, 2); % array name
	end
	m.p.stage = s; % store: 1 x ss, structure

	% Define plastic synapses
	m.p.synHebb = ["gen", "betaEx"; "gen", "betaIn"]; % Hebbian synapses
	m.p.synIn = ["betaIn", "betaEx"]; % syn. with spike-timing-independent plast.

	% Calculate the subcortical convergence function radii
	e = m.p.ecc; % current eccentricity (deg)
	m.p.radOpt = polyval(m.p.radOptCoef, e); % point spread function (deg)
	%	m.p.radHor = exp(polyval(m.p.radHorCoef, e)); % hor. cell rec. field (deg)
	m.p.radHor = model.readTab("radHor", e, m); % hor. cell rec. field (deg)
	m.p.radGang = model.readTab("radGang", e, m); % ganglion cell dendritic (deg)
	m.p.radCen = sqrt(m.p.radOpt ^ 2 + m.p.radGang ^ 2); % centre radius (deg)
	m.p.radSur = sqrt(m.p.radOpt ^ 2 + 2 * m.p.radHor ^ 2); % surround rad. (deg)

	% Calculate the cortical convergence function radii
	fac = model.readTab("ratCort", e); % inverse cortical mag'n factor (deg/mm)
	m.p.ratCort = fac; % store
	m.p.radBeta = fac * m.p.radBetaMm; % 4CBeta radius (deg)
	m.p.radInEx = fac * m.p.radInExMm; % inhibitory-excitatory radius (deg)
	m.p.radBetaCen = sqrt(m.p.radCen ^ 2 + m.p.radBeta ^ 2); % beta centre (deg)
	m.p.radBetaSur = sqrt(m.p.radSur ^ 2 + m.p.radBeta ^ 2); % beta surround (deg)

	% Calculate the arrays
	m = m.cone(m); % cone array
	m = m.gang(m); % ganglion cell array
	m = prim.gen(m); % geniculate array
	m = m.betaEx(m); % 4CBeta excitatory array
	m = m.betaIn(m); % 4CBeta inhibitory array

	%	Calculate the convergence functions
	m = prim.weight(m); % convergence functions

function s = stageList % create a list of stages and their arrays

	s = [
    "coneInt"    "cone"   
    "cone"       "cone"   
    "hor"        "cone"   
    "back"       "cone"   
    "biOff"      "cone"   
    "biOn"       "cone"   
    "gangOff"    "gangOff"
    "genOff"     "gangOff"
    "gangOn"     "gangOn" 
    "genOn"      "gangOn" 
		"gen"        "gen"    
    "betaIn"     "betaIn" 
    "sumGen"     "betaEx" 
    "sumIn"      "betaEx" 
    "betaEx"     "betaEx" 
		]; % ss x 2, string
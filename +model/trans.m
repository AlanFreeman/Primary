function p = trans(ec, pPriv, pEmp, m)

%	Purpose:
%		Manage transition of parameter (e.g. density or radius), p,
%		from private line to convergence
%	Input:
%		ec = required eccentricities (deg)
%		pPriv = private-line values
%		pEmp = empirical values
%		m = metadata
%	Ouput:
%		p = adjusted parameter

	%	Initialise
	eS = m.p.ecTrans(1); eE = m.p.ecTrans(2);
		% eccentricities at start and end of transition interval (deg)
	p = zeros(size(ec));
		
	%	Set parameter at private line eccentricities
	i = ec <= eS; % indices of private-line eccentricities
	p(i) = pPriv(i); % set parameter

	%	Set parameter in transition interval
	i = ec > eS & ec <= eE; % indices of eccentricities in transition interval
	f = (ec(i) - eS) ./ (eE - eS); % fraction of interval from start
	p(i) = (1 - f) .* pPriv(i) + f .* pEmp(i); % linear combination

	%	Set parameter at eccentricities beyond transition
	i = ec > eE; % indices of eccentricities beyond transition
	p(i) = pEmp(i); % set parameter
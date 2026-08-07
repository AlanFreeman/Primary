function d = scatFunc(d, m) % prepare scatter plot

	%	Obtain plot variable names
	n = string(m.scat.name); % variable names
	nX = n(1); nY = n(2); % names of horizontal, vertical coordinates
	nC = n(3); % name of colour variable

	%	Obtain plot variable values
	l = shiftdim(d.loc, 1); % cell location (deg): ls x 2
	x = d.resp.(nX).val; % x value: ls x 1
	y = d.resp.(nY).val; % y value: ls x 1

	%	Store and sort
	if nC ~= "" % colour available
		c = d.resp.(nC).val; % colour: ls x 1
		d = table(l, x, y, c, variableNames = ["loc", nX, nY, nC]);
		if isfield(m.scat, "descend") % change sort direction
			dir = 'descend'; % descending
		else
			dir = 'ascend'; % default
		end
		d = sortrows(d, nC, dir); % sort rows by colour
	else % no colour on which to sort
		c = repmat(missing, size(l, 1)); % no colour: ls x 1
		d = table(l, x, y, c, variableNames = ["loc", nX, nY, "col"]);
	end
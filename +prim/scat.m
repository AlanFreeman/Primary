function [d, m] = scat(d, m)

%	Prepare for scatter plot of neuronal indices
%
%	Inputs:
%		d.loc = neuronal locations (deg): 1 x ls x 2
%		d.resp.(n).val = value of index in m.scat.name, where n is the index name

	%	Obtain neuron locations and index names
	loc = shiftdim(d.loc, 1); % neuron location (deg): ls x 2
	n = string(m.scat.name); % index names: 1 x ns, string
	ns = length(n); % number of names

	%	Read index values
	v = cell(1, ns); % allocate storage
	for i = 1: ns % loop over indices
		nC = n(i); % name of current index
		v{i} = d.resp.(nC).val; % assign value: ls x vs
	end

	%	Store and sort
	d = table(loc, v{:}, variableNames = ["loc", n]); % store in new table
	if isfield(m.scat, "sort") % name of sort variable
		nS = m.scat.sort; % specified by user
	else
		nS = n(end); % default
	end
	if isfield(m.scat, "dir") % change sort direction
		dir = m.scat.dir; % specified by user
	else
		dir = "ascend"; % default
	end
	d = sortrows(d, nS, dir); % sort rows

	%	Choose neurons
	if isfield(m.scat, 'choose')
		x = d.loc(:, 1); y = d.loc(:, 2); % horizontal and vert. loc. (deg): ls x 1
		switch m.scat.choose
			case "cone" % weak spatial opponency, strong cone opponency
				i = d.soi < .05 & d.coiR > .85 & abs(x) < .05 & abs(y) < .05;
			case "def" % all indices are defined
				i = false(size(d, 1), 1); % initialise
				for j = 1: ns % loop over indices
					nC = n(j); % name of current index
					i = i | isnan(d.(nC)); % undefined indices
				end
				i = ~ i; % indices of defined indices
			case "spatial" %	strong spatial opponency, weak cone opponency
				i = d.soi > .85 & d.coiR < .25 & abs(x) < .05 & abs(y) < .05 & ...
					d.osiE > .4;
		end
		d = d(i, :); % select rows
	end
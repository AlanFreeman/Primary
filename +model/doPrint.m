function [d, m] = doPrint(d, m) % print variable for storage in readTab

	%	Initialise
	v = d.(m.print.var); % variable to print
	vs = length(v); % number of values to print
	if isfield(m.print, "vals") % number of values per line
		n = m.print.vals; % user specified
	else
		n = 8; % default
	end
	ls = ceil(vs / n); % number of lines to print

	%	Print
	for i = 1: ls % print a line
		j = n * (i - 1) + 1; % first index to print
		k = min(n * i, vs); % last index to print
		fprintf("%6.4f, ", v(j: k)); % print line
		fprintf("\n"); % finish line
	end
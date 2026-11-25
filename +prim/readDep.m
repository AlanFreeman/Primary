function [val, dim] = readDep(var, name) % read value of dependent variable

%	Input:
%		var = container for dependent variable: numerical or structure
%		name: required for structure:
%			target variable, for example, "dir.pref.val"
%	Output:
%		val = value of dependent variable
%		dim = dimensions of dependent variable:
%			it is assumed that they are in name, with "val" replaced by "dim",
%			for example "dir.pref.dim"

	val = var; % default
	if isstruct(var) % container is a structure

		%	Obtain value
		n = split(string(name), "."); % fields in name: 1 x ns, string
		ns = length(n); % number of fields
		for i = 1: ns % loop over fields
			val = val.(n(i)); % add fields recursively
		end

		%	Obtain dimensions
		dim = var; % initial name
		n(end) = "dim"; % replace "val by "dim"
		for i = 1: ns % loop over fields
			dim = dim.(n(i)); % add fields recursively
		end
		
	end
function d = readFolder(folder) % read raw data in folder and return table

	% Read contents file
	d = readtable([folder, '/', 'Contents'], 'textType', 'string', ...
		'delimiter', 'tab', 'consecutiveDelimitersRule', 'join'); % contents table
	f = d.file; % names of files
	c = d.content; % names of input variables
	l = d.label; % labels of output data
	v = split(d.var(1), ', ', 2); % names of output variables
	desc = split(d.desc(1), ', ', 2); % descriptions of output variables
	
	% Complile a table from each mat-file
	fs = length(f); % number of files
	dC = cell(fs, 1); % allocate storage
	for i = 1: fs % loop over files
		fC = f(i) + ".mat"; % name of current file
		s = load([folder, filesep] + fC); % file contents: structure
		cC = c(i); % current content
		cC = s.(cC); % file contents: double
		cs = size(cC, 1); % number of rows
		lC = repmat(l(i), [cs, 1]); % contents label
		dC{i} = table(lC, cC(:, 1), cC(:, 2), 'variableNames', ["label", v]);
	end

	% Combine and store
	d = vertcat(dC{:}); % combine data tables
	d.Properties.VariableDescriptions = ["label", desc]; % describe variables

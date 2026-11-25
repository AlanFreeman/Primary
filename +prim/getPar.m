function par = getPar(d, m, task) % get stimulus parameters and their values

%	Output:
%		par = structure with one member for each stimulus parameter:
%			par.name = parameter name: string
%			par.val = parameter value
%			par.vals = number of parameter values: double

	%	Initialise
	import prim.getIm; % find function
	if nargin > 2 % task is specified: get parameters from task metadata
		name = fieldnames(m.(task)); % names in m.task
	else % no task: get parameters from data table
		task = []; % no task
		name = d.Properties.CustomProperties.RespDim; % response parameters
	end
	
	% Find the stimulus parameters to be varied
	nameAll = fieldnames(m.p); % names of model parameters
	i = ismember(name, nameAll); % m.task fields that are stimulus par.
	name = string(name(i)); % names of specified parameters: 1 x ns, string
	pars = length(name); % number of parameters
	
	% Find the values of the stimulus parameters, which may be multi-column
	par(pars) = struct; % allocate storage
	for i = 1: pars % loop over stimulus parameters
		nameC = name(i); % name of current parameter
		par(i).name = nameC; % set name
		if ~ isempty(task) % task present
			valC = m.(task).(nameC); % parameter values
			par(i).val = valC; % set values
			par(i).vals = size(valC, 1); % set number of values
		else % no task
			valC = d.(nameC); % parameter values
			par(i).val = valC; % set values
			par(i).vals = size(valC, 2); % set number of values
		end
	end

	%	If the stimulus is an image, load the image file
	if m.p.stimS == "image" % *** fix ***
		%	[name, val, vals] = getIm(m, task); % get images to be used as stimuli
	end

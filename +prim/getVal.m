function [name, val, vals] = getVal(m, task)

% Get values of the stimulus parameters
%
%	Outputs:
%		name = names of specified stimulus parameters: cell
%		val = values for the specified stimulus parameters:
%			cell array with one member for each stimulus parameter
%		vals = number of values for each specified stimulus parameter: double

	switch m.p.stimS % grating or image
		case "grating"

			% Find the stimulus parameters to be varied
			name = fieldnames(m.(task)); % names in m.task
			def = prim.setDef; % default parameters
			nameAll = fieldnames(def.p); % names of default parameters
			i = ismember(name, nameAll); % m.task fields that are stimulus parameters
			name = name(i)'; % names of specified parameters: 1 x ns
			names = length(name); % number of parameters
			
			% Find the values of the response parameters, which may be multi-column
			val = cell(1, names); % values of parameters
			vals = zeros(1, names); % number of values
			for i = 1: names % loop over response parameters
				nameC = name{i}; % name of current parameter
				valC = m.(task).(nameC); % values of parameter
				val{i} = valC; % store values
				vals(i) = size(valC, 1); % update number of values
			end

		case "image"
			[name, val, vals] = prim.getIm(m, task); % get images to be used as stim.
	end
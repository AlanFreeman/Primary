function [d, m] = stage(d, m) % list properties of specified stages

%	Assumptions:
%		input table is one row of the modulation file *** is this necessary, useful?
%		m.stage.name = stage names: cell

	%	Construct structure of stages
	name = string(m.stage.name); % names of stages: string
	name = name(:); % make it a column vector: ss x 1
	ss = length(name); % number of stages
	p = cell(ss, 1); % empty array: ss x 1, cell
	p = struct('name', p); % stage properties: ss x 1, struct
	for i = 1: ss % loop over stages

		%	Assemble data for current stage
		nameC = name(i); % name of current stage
		j = ismember([m.p.stage.name], nameC); % index in list of stages
		s = m.p.stage(j); % structure for current stage
		a = s.array; % name of array to which stage belongs

		%	Store properties for current stage
		p(i).name = nameC; % stage name: string
		p(i).loc = m.p.(a).loc; % cell locations (deg): ls x 2
		p(i).type = m.p.(a).type; % cell type

	end

	%	Store
	d = repmat(d, [ss, 1]); % one row for each stage
	d.stage = name; % stage name
	d.prop = p; % stage properties
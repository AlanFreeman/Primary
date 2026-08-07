function out = setFig(m, opt) % set metadata for figure plotting

	arguments

		m struct % metadata
		opt.fig1
		opt.fig2
		opt.figIndexB
		opt.init
		opt.xSym (1, 1) double
		opt.ySym (1, 1) double
		opt.zSym (1, 1) double
		opt.xZero (1, 1) double
		opt.yZero (1, 1) double
		opt.zZero (1, 1) double
		opt.zoom (1, 2) cell
		%	opt.lineWidth (1,1) {mustBeNumeric} = 1

	end

	%	Loop over options
	name = string(fieldnames(opt)); % option field names
	names = length(name); % number of options
	for i = 1: names
		nameC = name(i); % name of current option
		switch nameC
			case "fig1" % SfN 25 poster: frequencies for spatial frequency response
				freqS = linspace(.3, 30, 10);
				out = freqS';
			case "fig2" % SfN 25 poster: stimuli for stimulus sequence
				switch opt.(nameC) % stimulus number
					case 1, m.p.dir = 0; m.p.freqS = 8;
					case 2, m.p.dir = 60; m.p.freqS = 12;
					case 3, m.p.dir = -30; m.p.freqS = 4;
				end
				out = m;
			case "figIndexB" % Prim paper
				lim = opt.(nameC); % maximum amplitude (Hz)
				out{1} = [0, lim]; % for clim
				out{2} = {'ticks', lim * [0, .5, 1]}; % for colorbar
			case "init"

				%	Default variables
				m.x = "x"; % default independent variable
				m.y = "resp"; % default independent or dependent variable
				m.z = "resp"; % default dependent variable

				%	Default setting for plotting variables
				m.pub.units = "centimeters"; % define unit before using it
				m.pub.axisLength = [5, 5]; % length of x-, y-axes (cm)
				m.pub.export = {"colorSpace", "cmyk"}; % for file export
				m.pub.lineWidth = 1.5; % width of axes and lines (cm)
				out = m;

			case {"xSym", "ySym", "zSym"} % limits and ticks symmetric about zero
				o = extract(nameC, 1); % first character in name
				val = opt.(nameC); % option value
				limSym = .5 * [-1, 1]; tickSym = .5 * [-1, 0, 1];
					% limits and ticks symmetric about zero
				out{i} = {o + "Lim", val * limSym, o + "Tick", val * tickSym};
					% store
			case {"xZero", "yZero", "zZero"} % limits and ticks starting at zero
				o = extract(nameC, 1); % first character in name
				val = opt.(nameC); % axis span (deg)
				limZero = [0, 1]; tickZero = [0, .5, 1];
					% limits and ticks starting at zero
				out{i} = {o + "Lim", val * limZero, o + "Tick", val * tickZero}; % out
			case "zoom" % zoom into cell of interest
				val = opt.(nameC); % option value: 1 x 2, cell
				loc = val{1}; % cell location (deg)
				span = val{2}; % reduction in visual field width
				span = span * .5 * m.p.wid * [-1, 1];
				xLim = loc(1) + span; yLim = loc(2) + span;
				out{i} = {"xLim", xLim, "yLim", yLim};
		end
	end

	%	Store output
	if iscell(out)
		out = [out{:}]; % combine strings
	end

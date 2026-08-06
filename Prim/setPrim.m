function setPrim % calculate Prim model constants from published work

	% Initialise
	import prim.readTab prim.setLit; % function location
	m = setPrimTask; % define task function handles
	m = setLit(m); % set metadata calculated from the literature
	m.project = 'Prim'; % results are stored in userpath/Data/m.project
	m.folder = [userpath, '/Data/Model']; % data folder

	% Determine and display parameters
	switch 'dens.ecc gang' % select analysis tasks
		case 'conv' % convert from XYZ colour space to LMS colour space
			m.tasks = 'xyz2lms';
		case 'count.ecc cone' % cone count vs eccentricity: Packer (89)
			m.tasks = 'densCone countCone plot set add';
			m.x = 'eccDeg'; m.y = 'count';
			m.set.axes = {'xLim', [0, 60]};
			m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
				'xTick', [.1, 1, 10, 60], 'yScale', 'log'};
			m.add.funFun = @(d, m) plot(1, interp1(d.eccDeg, d.count, 1), 'ok');
		case 'count.ecc cen' % % number of cones per ganglion cell centre
			m.tasks = 'densCone select countCen plot set';
			m.select.quad = 'temporal';
			m.x = 'eccDeg'; m.y = 'count';
			m.set.axes = {'xScale', 'log', 'xLim', [.01, 40], 'xTick', ...
				[.01, .1, 1, 10], 'yLim', [0, 20], 'yTick', [0, 10, 20], ...
				'clipping', 'on'};
		case 'count.ecc dens' % cones per ganglion cell via functional g.c. density
			m.tasks = 'densGang countGang densGangFun densGangSub ratDens';
			m.tasks = [m.tasks, ' ', 'plot set add'];
			m.x = 'ecc'; m.y = 'ratio';
			m.ratDens.group = 'type';
			m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
				'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [.3, 20], ...
				'yTick', [.5, 1, 10]};
			m.add.funFun = @(d, m) plot(.5, 4484.03 / 9275.5, 'ok');
				% from ganglion cell count
		case 'count.ecc gang' % ganglion cell count vs functional ecc.: Wässle (89)
			m.tasks = 'densGang countGang plot set add';
			m.y = 'count';
			switch 'func' % type of eccentricity
				case 'anat' % anatomical
					m.tasks = 'densGang countGang plot set add';
					m.x = 'eccAnat';
					m.set.axes = {'xScale', 'log', 'xLim', [1, 60], ...
						'yScale', 'log', 'yLim', [100, 1e7], 'clipping', 'off'};
					e = m.p.magRet * .55 / .9; n = 33000; % data for 550 um circle
					m.add.funFun = @ (d, m) plot(e, n, 'o');
				case 'func' % functional
					m.x = 'ecc';
					switch 'log' % linear or log?
						case 'lin' % linear
							lim = 1;
							m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
								'yLim', [0, 4e4], 'yTick', 4e4 * [0, .5, 1]};
						case 'log' % logarithmic
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 100], ...
								'yScale', 'log', 'yLim', [100, 1e7], 'clipping', 'off'};
					end
			end
			m.add.funFun = @(d, m) plot(1, interp1(d.ecc, d.count, 1), 'ok');
		case 'count.wave' % histogram of preferred wavelengths in natural images
			m.tasks = 'hist plot set';
			m.folder = [userpath, '/Data/Prim']; m.p.file = 'Image 10';
			m.x = 'cen'; m.y = 'count';
			m.hist.group = 'name'; m.hist.x = 'wave.wave';
			m.plot.group = 'name';
		case 'cover.ecc' % ganglion cell coverage
			m.tasks = ['densGang countGang densGangFun densGangSub ', ...
				'select cover plot set'];
			m.x = 'ecc'; m.y = 'cover';
			m.select.type = 'mid';
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yLim', [0, 5], 'yTick', [0, 2.5, 5]};
			m.set.axes = {'xScale', 'log', 'xLim', [.01, 40], ...
				'xTick', [.01, .1, 1, 10, 40], 'yLim', [0, 15], 'yTick', [0, 7.5, 15]};
		case 'dens beta' % 4CBeta cell density: O'Kusky (82), Fitzpatrick (87)
			switch 'area' % method
				case 'area' % areal densities
					a = [330, 445, 516, 614, 175, 314, 381, 393]; % all neurons: Fitz (87)
					i = [43, 86, 79, 92, 31, 67, 55, 52]; % GABA neurons
					r = sum(i) / sum(a); % ratio of inhib. to all neurons in 2, 3, 4CBeta
					dens = 30e3; % density of 4CBeta (cells/mm^2): O'Kusky (82)
					unit = "mm^-2"; % spatial unit
				case 'volume' % volumetric densities
					r = .155; % ratio of inhib. to all neurons in 4CBeta: Fitz (87)
					dens = 211100; % cell density in 4CBeta, cells/mm^3: O'Kusky (82)
					unit = "mm^-3"; % spatial unit
			end
			m.p.densBeta = [1 - r, r] * dens; % density of ex., inhibitory cells
			fprintf("Density of exc., inh. cells in 4CBeta (" + unit + "): " + ...
				"%g, %g\n", m.p.densBeta); % display
			return			
		case 'dens.ecc cone' % cone density: Packer (89)
			m.tasks = 'densCone select prep plot set add';
			field = 'both'; % visual field: 'both', 'far', 'near'
			save = 0; % save result: 0 or 1
			m.x = 'eccDeg'; m.y = 'densDeg';
			m.select.quad = 'temporal';
			m.plot.group = 'quad'; m.plot.arg = {'-'};
			m.add.funFun = @(d, m) plot(.5, 14087 / pi, 'ok'); % from cone count
			switch field % plotting commands
				case 'both' % both far and near
					yLim = [3, 2e4];
					switch 'lin' % linear or logarithmic x axes
						case 'lin' % linear
							xLim = 55; xTick = 50;
							m.set.axes = {'xLim', xLim * [0, 1], 'xTick', ...
								xTick * [0, .5, 1], 'yScale', 'log', 'yLim', yLim};
						case 'log' % log
							z = .3; m.prep.zero = z; % foveal proxy
							m.set.axes = {'xScale', 'log', 'xLim', [z, 60], ...
								'xTick', [z, 1, 10, 60], 'yScale', 'log', 'yLim', yLim, ...
								'clipping', 'on'};
					end
				case 'far' % as published
					m.select.near = 0;
					m.x = 'ecc'; m.y = 'dens'; m.plot.arg = {'-s'}; lim = 1000;
					m.set.axes = {'xLim', [0, 20], 'xTick', [0, 10, 20], ...
						'yLim', lim * [1, 25], 'yTick', lim * [5, 15, 25]};
				case 'near' % as published
					m.select.near = 1;
					m.x = 'ecc'; m.y = 'dens'; m.plot.arg = {'-s'}; lim = 1000;
					m.set.axes = {'xLim', [0, 1.25], 'xTick', [.25, .75, 1.25], ...
						'yLim', lim * [25, 225], 'yTick', lim * [25, 125, 225]};
			end
			switch save % save results? *** check ***
				case 1
					m.tasks = 'densCone select doPrint'; 
					m.doPrint.var = 'densDeg';
			end
		case 'dens.ecc gang' % ganglion cell density: Wässle (89)
			m.tasks = 'densGang select plot set';
			m.select.type = {'anatNarr', 'anatWide'};
			switch 'deg' % which units?
				case 'deg' % convert to deg
					m.x = 'eccDeg'; m.y = 'densDeg';
					yLim = [3, 2e4];
					switch 'log' % linear or logarithmic eccentricity?
						case 'lin' % linear
							m.set.axes = {'xLim', [0, 55], 'xTick', [0, 25, 50], ...
								'yScale', 'log', 'yLim', yLim};
						case 'log' % logarithmic
							m.set.axes = {'xLim', [.3, 60], 'xTick', [.3, 1, 10, 60], ...
								'xScale', 'log', 'yScale', 'log', 'yLim', yLim, ...
								'clipping', 'off'};
					end
				case 'mm' % plot as published
					m.x = 'ecc'; m.y = 'dens';
					xLim = 12; xTick = [0, 6, 12]; yLim = [100, 1e6];
					m.set.axes = {'xLim', xLim * [0, 1], 'xTick', xTick, ...
						'yScale', 'log', 'yLim', yLim};
			end
			m.plot.line = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
		case 'dens.ecc gang fit' % fit ganglion cell density: Wässle (89)
			switch 'show'
				case 'check'
					m.tasks = 'densGang select plot set add';
					m.add.funFun = @(d, m)plot(d.eccDegLog, ...
						exp(polyval(m.p.densGangCoef, d.eccDegLog))); % *** use readTab ***
				case 'fit'
					m.tasks = 'densGang select plot set fitglm pred add';
				case 'plot'
					m.tasks = 'densGang select plot set';
				case 'show'
					m.tasks = 'densGang select fitglm show'; % set m.p.densGangCoef here
			end
			m.select.type = {'anatNarr', 'anatWide'};
			m.x = 'eccDegLog'; m.y = 'densDeg';
			m.plot.line = 'type';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', log10([.3, 60]), 'yScale', 'log', ...
				'yLim', [3, 2e4], 'clipping', 'off'};
			m.fitglm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y, ...
						'link', 'log'};
			m.add.line = {};
		case 'dens.ecc gang fitlm' % fit ganglion cell density: Wässle (89)
			m.tasks = 'densGang select plot set fitlm show pred add';
			m.select.type = {'anatNarr', 'anatWide'};
			m.x = 'eccDegLog'; m.y = 'densDegLog';
			%	m.plot.line = 'type';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
				'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4], ...
				'clipping', 'off'}; % *** fix ***
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};			
		case 'dens.ecc gang fun' % g.c. density with functional ecc.: Wässle (89)
			m.tasks = 'densGang countGang densGangFun prep plot set';
			m.x = 'ecc'; m.y = 'densDeg';
			switch 'log' % x axis scaling
				case 'lin' % linear
					xLim = 55; yLim = [3, 3e4];
					m.set.axes = {'xLim', xLim * [0, 1], 'xTick', [0, 25, 50], ...
						'yScale', 'log', 'yLim', yLim, 'clipping', 'off'};
				case 'log' % logarithmic
					m.prep.zero = .1;
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
						'xTick', [.1, 1, 10, 60], 'yScale', 'log', 'yLim', [4, 2e4], ...
						'clipping', 'off'};
			end
		case 'dens.ecc gang sub' % g.c. dens. subpopulations: multiple authors
			m.tasks = 'densGang countGang densGangFun densGangSub prep plot set';
			%	m.tasks = ...
			%		'densGang countGang densGangFun densGangSub prep plot set add';
			%	m.tasks = ['densGang countGang densGangFun densGangSub ' ...
			%	'select priv doPrint']; % use this for readTab
			z = .3; m.prep.zero = z;
			m.select.type = 'mid';
			m.x = 'ecc'; m.y = 'densDeg';
			m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [z, 60], 'xTick', ...
				[z, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4]};
			m.add.funFun = @(d, m) plot(.5, (34898 / pi) * [1, m.p.ratGang], ...
				'o'); % from ganglion cell count
			m.priv.var = 'densGang'; m.doPrint.var = 'densOn';
		case 'eccFun.ecc' % func. ecc. vs anat. ecc.: McGregor (18), Schein (88)
			m.tasks = 'eccFun plot set add'; % plot
			m.x = 'ecc'; m.y = 'eccFun';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			lim = 15;
			m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
				'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'on'};
			m.add.funFun = @(d, m)plot([0, 90], [0, 90], '--');
		case 'eccFun.ecc fit' % fun. vs anat. ecc. + fit: McGregor (18), Schein (88)
			m.tasks = 'eccFun plot set eccFunFit add'; % plot
			%	m.tasks = 'eccFun eccFunFit save'; % save
			%	m.eccFunFit.mcGregor = .03; % lower limit on McGregor fun. ecc.
			m.eccFunFit.schein = 5; % lower limit on Schein eccentricities
			m.x = 'ecc'; m.y = 'eccFun';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			switch 'log' % axis scaling
				case 'lin' % linear
					lim = 10;
					m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
						'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'on'};
				case 'log' % logarithmic
					m.set.axes = ...
						{'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
						'yScale', 'log', 'yLim', [.007, 10], 'yTick', [.01, .1, 1, 10], ...
						'clipping', 'on'};
			end
			m.add.funFun = @(d, m) plot(d.ecc, d.eccFun, 'k');
			m.save.name = 'Functional eccentricity';
		case 'eccFun.ecc mcG' % func. ecc. vs anat. ecc.: McGregor (18)
			m.tasks = 'eccFunMcG plot set';
			m.x = 'eccFunUm'; m.y = 'eccUm'; m.plot.arg = {'o'};
			m.set.axes = {'xLim', [0, 140], 'xTick', [0, 70, 140], ...
				'yLim', [100, 600], 'yTick', [100, 300, 500]};
		case 'freq.ecc' % calc. spat. freq. cutoff of 4CBeta cells: pragmatic
			m.tasks = 'freq'; % print results and enter in readTab
			m.x = 'ecc'; m.y = 'freq';
			m.freq.cutoff = .1; % attenuation of spat. freq. resp. at cutoff
		case 'gain' % surround gain: Croner (95)
			freq = 2 * pi * 4.22; % stimulus temporal frequency (radians/s)
			a = 1 + 1i * m.p.tau * freq; % temporal attenuation
			r = .547; % surround response / centre response at 0 spatial frequency
			switch 'fzero' % method
				case 'direct' % algebraic solution
					k = (a ^ 2) * ((1 / r) - 1);
					fprintf('Surround gain magnitude, phase (deg): %3g, %3g\n', ...
						abs(k), (180 / pi) * angle(k)); % list
				case 'fzero' % numerical solution
					fun = @ (k) abs(1 + k / a ^ 2) - 1 / (1 - r); % find zero for this fun
					m.p.kSur = fzero(fun, 1); % solve for kSur
					fprintf('Surround gain: %3g\n', m.p.kSur); % list
			end
			return
		case 'image' % process images from ImageNet: Deng (09)
			m.tasks = 'procIm save';
			m.procIm.files = 1001: 1591; % serial numbers of image files
			m.procIm.proc = 'filter'; % type of processing: filter or screen
			%	m.save.folder = [userpath, '/Data/Prim'];
			m.save.name = 'Image 7.mat';
		case 'image.x.y' % show images
			m.tasks = 'showIm plot set';
			m.x = 'x'; m.y = 'y'; m.z = 'z';
			m.showIm.files = 367;
			m.plot.group = 'name';
			m.plot.funFun = @(d, m) image(d.x, d.y, flipud(d.image{1}));
			m.set.axes = {'xLim', [0, 500], 'yLim', [0, 500], 'clipping', 'off'};
			m.set.colorbar = {};
			m.set.title = {'Interpreter', 'none'};
		case 'mtf.freq' % modulation transfer function: Navarro (93)
			m.tasks = 'mtf plot set';
			m.x = 'freq'; m.y = 'mtf';
			m.plot.line = 'ecc';
			m.set.axes = {'yLim', [.005, 1], 'yScale', 'log'};			
		case 'offset.ecc' % ganglion cell offset vs. eccentricity: Schein (88)
			m.tasks = 'offset plot set';
			m.plot.line = 'source'; m.plot.arg = {'clipping', 'off'};
			m.x = 'eccMm'; m.y = 'offsetMm'; % m.plot.arg = {'o'};
			m.set.axes = {'xLim', [0, 2.5], 'xTick', 0: 2, ...
				'yLim', [0, .4], 'yTick', [0, .2, .4]};
		case 'psf.loc' % point spread function: Navarro (93)
			m.tasks = 'mtf psf plot set';
			switch 'psf' % function to display
				case 'mtf' % modulation transfer function
					m.x = 'freq'; m.y = 'mtf';
					m.plot.line = 'ecc';
					m.set.axes = {'xLim', [0, 60], 'yScale', 'log', 'yLim', [.005, 1]};
				case 'psf' % point spread function
					m.x = 'loc'; m.y = 'psf';
					m.plot.line = 'ecc';
			end
		case 'rad beta' % 4CBeta radius: multiple authors
			m.tasks = 'radBeta plot set'; % individual radii
			m.tasks = 'radBeta plot set list'; % m.p.radBetaMm = convergence radius
			m.x = 'loc'; m.y = 'radius'; m.plot.arg = {'o'};
			m.plot.line = 'source';
			m.set.line = 'fill';
		case 'rad.ecc' % centre, surround radius decomposition: multiple authors
			m.x = 'ecc'; m.y = 'radius';
			m.plot.line = 'source';
			m.add.arg = {'o'};
			x = {'xScale', 'log', 'xLim', [.006, 40], 'xTick', [.01, .1, 1, 10]};
			switch 'cen' % choose mechanism
				case 'cen' % centre mechanism
					switch 'mech' % which display?
						case 'emp' % empirical data
							m.tasks = 'rad select plot set radCen add'; % empirical data
							m.select.source = {'cenConv'};
						case 'mech' % mechanisms
							m.tasks = 'rad select plot set'; % mechanisms
							m.select.source = {'optics', 'denGang', 'cenConv'};
					end
					y = {'yScale', 'log', 'yLim', [.007, .3], 'yTick', [.01, .1]};
					m.set.axes = [x, y];
				case 'sur' % surround mechanism
					switch 'emp' % choose plot
						case 'emp' % empirical data
							m.tasks = 'rad select plot set radSur add'; % empirical data
							m.select.source = {'surConv'};
						case 'mech' % mechanisms
							m.tasks = 'rad select plot set'; % mechanisms
							m.select.source = {'optics', 'fieldHor', 'surConv'};
					end
					yLim = [.007, 5]; yTick = [.01, .1, 1];
					m.set.axes = {'xScale', 'log', 'xLim', [.006, 40], ...
						'xTick', [.01, .1, 1, 10], 'yScale', 'log', 'yLim', yLim, ...
						'yTick', yTick, 'clipping', 'off'};
			end
		case 'rad.ecc beta' % geniculocortical convergence radius: multiple authors
			m.tasks = 'radBetaEcc plot set';
			m.p.file = 'Cortical magnification'; % *** use readTab? ***
			m.x = 'ecc'; m.y = 'radius';
			m.set.axes = {'xLim', [0, 50], 'xTick', [0, 25, 50]};
		case 'rad.ecc cen' % gang. cell centre radius: Croner (95), Lee (98)
			m.tasks = 'radCen select plot set'; % no fitting
			m.select.source = 'Lee';
			m.x = 'ecc'; m.y = 'radius';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			switch m.select.source % plot in published form
				case 'Croner' % Croner (95)
					m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40], 'yLim', [0, .3]};
				case 'Lee' % Lee (98)
					m.y = 'dev';
					m.set.axes = {'xLim', [0, 15], 'xTick', [0, 7.5, 15], ...
						'yScale', 'log', 'yLim', [.1, 10]};
			end
		case 'rad.ecc cen fit' % fit gang. cell centre radius vs ecc.: Croner (95)
			m.tasks = 'radCen select plot fitlm pred add set'; % show prediction
			m.tasks = 'radCen select fitlm show'; % show model
			m.tasks = 'radCen select plot add set'; % check stored parameters
			m.select.source = 'Croner';
			m.x = 'ecc'; m.y = 'radiusLog';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', 'eccFun', ...
				'responseVar', 'radiusLog'};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.radCen; m.add.funFun = @(d, m)plot(d.ecc, p(1) + ...
					p(2) * d.ecc + p(3) * d.ecc .^ 2 + p(4) * d.ecc .^ 3);
			end
			m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40]};
		case 'rad.ecc cone' % cone inner segment diameter: Packer (89)
			m.tasks = 'radCone';
		case 'rad.ecc gang' % ganglion cell dendritic diameter: Watanabe (89)
			m.tasks = 'radGang plot set'; % plot as published
			m.x = 'eccMm'; m.y = 'diam';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', [0, 14], 'xTick', [0, 7, 14], ...
				'yScale', 'log', 'yLim', [1, 500]};		
		case 'rad.ecc gang field' % gang. cell cen., sur. radius: Godat (22)
			m.tasks = 'radGangField plot set';
			m.x = 'ecc'; m.y = 'radius';
			m.plot.arg = {'o'}; m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [.006, .1], 'xTick', [.01, .1], ...
				'yScale', 'log', 'yLim', [.01, .1]};
		case 'rad.ecc gang fit' % ganglion cell dendritic diameter: fitted curve
			m.x = 'eccFun'; m.y = 'radius';
			switch 'check' % select analysis tasks
				case 'check' % check stored parameters
					%	m.tasks = 'radGang plot set add'; % use polynomial
					%	m.add.funFun = @(d, m)plot(d.eccFun, ...
					%	polyval(m.p.radGangCoef, d.eccFun));
					m.tasks = 'radGang plot set fitlm priv add'; % use readTab
					m.add.funFun = @(d, m)plot(d.eccFun, readTab('radGang', d.eccFun));
					m.priv.var = 'radGang';
				case 'fit' % fit regression model
					m.tasks = 'radGang plot set fitlm pred add'; % add fitted curve
				case 'priv' % private line from cone to parafoveal ganglion cell
					m.tasks = 'radGang plot set fitlm priv add'; % add fitted curve
					m.tasks = 'radGang plot set fitlm priv doPrint'; % for readTab
					m.priv.var = 'radGang'; m.doPrint.var = 'radius';
				case 'show' % show regression model
					m.tasks = 'radGang fitlm show'; % model: set m.p.radGangCoef from this
			end
			m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};
			switch 'lin'
				case 'lin'
					m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
					'yLim', [.001, .4]};
					%	'yScale', 'log', 'yLim', [.001, .4]};
				case 'log'
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
						'xTick', [.1, 1, 10, 60], 'yScale', 'log', ...
						'yLim', [.007, .4], 'clipping', 'off'};
			end
		case 'rad.ecc hor dend' % hor. cell dendrite area: Wässle (89) Horizontal
			m.tasks = 'radHorDend plot set'; % plot as published
			m.x = 'ecc'; m.y = 'area';
			m.plot.group = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', [0, 14], 'xTick', [0, 7, 14], ...
				'yScale', 'log', 'yLim', [1e-4, 2e-2]};			
		case 'rad.ecc hor dend fit' % horizontal cell dendrite area: fitted curve
			m.tasks = 'radHorDend plot fitlm pred add set'; % add fitted curve
			%	m.tasks = 'radHorDend fitlm show'; % show regression model
			%	m.tasks = 'radHorDend plot add set'; % check stored parameters
			m.x = 'eccFun'; m.y = 'areaLog';
			m.plot.group = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
			m.fitlm.group = 'type'; m.pred.group = 'type';
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.area; m.add.funFun = @(d, m)plot(d.eccFun, p(1) + ...
					p(2) * d.eccFun + p(3) * d.eccFun .^ 2 + p(4) * d.eccFun .^ 3);
			end
			m.set.axes = {'xLim', [0, 70], 'xTick', [0, 35, 70], ...
				'yLim', [-4, log10(.02)]};			
		case 'rad.ecc hor field' % horizontal cell receptive field: Packer (02)
			m.tasks = 'radHorField select plot set'; % no fitting
			m.x = 'ecc'; m.y = 'diam';
			m.plot.arg = {'o'}; m.plot.line = 'source';
			m.set.axes = {'xLim', [0, 16], 'xTick', [0, 8, 16], ...
				'yLim', [0, 1200], 'yTick', [0, 400, 800]};
		case 'rad.ecc hor field fit' % hor. cell radius: Wässle (89), Packer (02)
			switch 'fit' % select analysis tasks
				case 'check' % check stored parameters
					m.tasks = 'radHor select plot set add';
					m.select.source = {'dend', 'conv'};
					m.add.funFun = @(d, m)plot(d.eccFun, ...
						exp(polyval(m.p.radHorCoef, d.eccFun)));
				case 'fit' % fit regression model
					m.tasks = 'radHor select plot set fitglm pred add';
					m.select.source = {'dend', 'conv'};
				case 'plot' % no fit
					m.tasks = 'radHor select plot set';
					m.select.source = {'dend', 'wide'};
				case 'show' % show regression model
					m.tasks = 'radHor select fitglm show'; % set m.p.radHorCoef from list
					m.select.source = {'dend', 'conv'};
			end
			m.x = 'eccFun'; m.y = 'radius';
			m.plot.arg = {'o'}; m.plot.line = 'source';
			m.fitglm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y, ...
				'link', 'log'};
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yScale', 'log', 'yLim', [.035, 1.3]};
		case 'rad.ecc opt' % point spread function: Navarro (93), Williams (81)
			m.tasks = 'radOpt plot set';
			switch 'check' % choose tasks to perform
				case 'add' % add fitted model to empirical data
					m.tasks = 'radOpt plot set fitlm pred add';
				case 'check' % check stored coefficients for fitted model
					m.tasks = 'radOpt plot set add';
					m.add.funFun = @(d, m)plot(d.ecc, polyval(m.p.radOptCoef, d.ecc));
				case 'show' % show statistics for fitted model
					m.tasks = 'radOpt plot set fitlm show'; % set m.p.radOptCoef from list
			end
			m.x = 'ecc'; m.y = 'radius'; m.plot.arg = {'-o'};
			m.fitlm.arg = {'radius ~ ecc^2 - ecc'};
			switch 'psf' % data to display
				case 'rrf' % retinal resolution function, as published
					m.y = 'rrf';
					m.set.axes = {'xLim', [-60, 60], 'xTick', [-60, 0, 60], ...
						'yScale', 'log', 'yLim', [.7, 100], 'yTick', [1, 10, 100]};
				case 'psf' % radius of point spread function
					m.set.axes = {'xLim', [-60, 60], 'xTick', [-60, 0, 60], ...
						'yLim', [0, .075], 'yTick', [0, .025, .05]};
			end
		case 'rad.ecc sur' % ganglion cell surround radius: Croner (95), Lee (98)
			m.tasks = 'radSur select plot set'; % no fitting
			m.x = 'eccFun'; m.y = 'radius';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
				'yScale', 'log', 'yLim', [.01, 10], 'clipping', 'off'};
		case 'rad.ecc sur fit' % ganglion cell surround radius: fitted curve
			% *** untested code ***
			m.tasks = 'radSur select plot fitlm pred add set'; % show prediction
			m.tasks = 'radSur select fitlm show'; % show model
			m.tasks = 'radSur select plot add set'; % check stored parameters
			m.select.source = 'sur';
			m.x = 'ecc'; m.y = 'radiusLog';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', 'ecc', ...
				'responseVar', 'radiusLog'};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.radSur; m.add.funFun = @(d, m)plot(d.ecc, p(1) + ...
					p(2) * d.ecc + p(3) * d.ecc .^ 2 + p(4) * d.ecc .^ 3);
			end
			m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40]};
		case 'rad fix' % radius of fixation eye movements: Skavenski (75)
			hor = [5.4, 4.6, 4.2, 4.8]; vert = [4.2, 6.4, 5.8, 7.8];
				%	standard deviations of the horizontal and vertical eye position (min)
			r = mean([hor, vert]) / 60; % mean standard deviation (deg)
			r = sqrt(2) * r; % radius (deg)
			fprintf('Radius of fixation eye movements, radFix = %4g\n', r); % report
			return
		case 'rat cone' % cone ratio: munds (22)
			switch 'munds' % choose source
				case 'mult' % multiple authors
					l = .6; % L / (L + M): Dacey et al (00) "Physiology ...", temporal q.
					lm = [l, 1 - l]; % [L, M] / (L + M)
					s(1) = .09934; % S / (L + M): Martin et al. (99); ecc. = 4-67 deg temp
					s(1) = s(1) / (1 + s(1)); % S / (L + M + S): Martin et al. (99)
					s(2) = .073; % S / (L + M + S): Roorda et al. (01); ecc. = 1-1.5 deg
					s = mean(s); % average
					r = [(1 - s) * lm(1), (1 - s) * lm(2), s]; % L: M: S ratio
				case 'munds' % Munds (22)
					lm = 1.03; % L / M
					s = .143; % S / (L + M)
					r = [lm / (1 + lm), 1 / (1 + lm), s] / (1 + s);
						% [L, M, S] / (L + M + S), m.p.ratCone = r
			end
			fprintf('Cone ratio, [L, M, S]: [%4g, %4g, %4g]\n', r); % report
			return
		case 'rat gang' % ratio of midget to all ganglion cells: multiple authors
			switch 'peng' % source
				case 'mult' % multiple authors: off-midget / all midget ganglion cells
					r = [.63, .62, .53]; % Dacey (93), Peng (19), Rhoades (19)
					r = mean(r); % average over studies
					fprintf('Off-midget / all midget ganglion cells: %4g\n', r); % report
				case 'peng' % Peng (19)
					r = [.47, .37; .5, .33]; % ratio of midget to all ganglion cells
						% row 1: fovea; row 2: periphery;
						% column 1: off-centre; column 2: on-centre
					r = mean(r); % mean of fovea, periphery
					fprintf('Ratio of off- to all midget g.c., m.p.ratSign: %g\n', ...
						r(1) / sum(r));
					fprintf('Ratio of midget to all g.c., m.p.ratGang: %g\n', sum(r));
			end
			return
		case 'rat gang dev' % nearest neighbour gang. cells: s.d./dist., Dacey (93)
			dist = 145; % mean distance between nearest neighbour on-centre g.c. (um)
			dev = 18; % standard deviation of distance (um)
			r = dev / dist; % ratio, m.p.kGangDev = r
			fprintf('Nearest neighbour gang. cells distance: s.d. / dist. = %g\n', r);
			return
		case 'rat ret' % retinal magnification ratio: Perry, Cowey (85)
			% Also see De Monasterio et al. (85): 1 / .203 = 4.93
			f(1) = 1 / .223; % M. mulatta
			f(2) = 1 / .201; % M. fascicularis
			f = mean(f); % average over species (deg/mm): set m.p.magRet = f
			fprintf('Retinal magnification factor (deg / mm): %4g\n', f); % report
			return
		case 'rat.ecc cort' % inverse cortical magnification factor: Dow (81)
			m.tasks = 'ratCort plot set';
			%	m.tasks = 'ratCort plot set save'; % don't use
			%	m.tasks = 'ratCort plot set doPrint'; % use this to save in readTab
			m.doPrint.var = 'fac';
			switch 'min'
				case 'deg' % plot in deg
					m.x = 'ecc'; m.y = 'fac';
					%	lim = 10;
					lim = 1; % for close view of foveal representation
					m.set.axes = {'xLim', lim * [0, 5], 'xTick', lim * [0, 2.5, 5] ...
						'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'off'};
					%	m.save.name = 'Cortical magnification';
				case 'min' % plot as published
					m.x = 'eccMin'; m.y = 'facMin';
					m.ratCort.eccs = 200;
					m.set.axes = {'xLim', [3, 3000], 'yLim', [3, 1000], ...
						'xScale', 'log', 'yScale', 'log', 'clipping', 'off'};
			end
		case 'rat.ecc cort fit' % inverse cortical magnification factor: Dow (81)
			% *** delete? ***
			switch 'spline' % regression method
				case 'log' % fit log factor versus log eccentricity
					m.tasks = 'ratCort plot set fitglm show'; % show regression model
					m.tasks = 'ratCort plot set fitglm pred add'; % fit regression model
					m.x = 'eccLog'; m.y = 'fac';
					m.plot.arg = {'o'};
					m.set.axes = {'yLim', [.05, 20], 'yScale', 'log', 'clipping', 'off'};
					m.fitglm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y, ...
						'link', 'log'};
				case 'spline' % interpolate with spline
					m.tasks = 'ratCort plot set stop pred add'; % fit regression model
					m.x = 'ecc'; m.y = 'fac';
					lim = 1; % lim = 10;
					m.set.axes = {'xLim', lim * [0, 5], 'xTick', lim * [0, 2.5, 5] ...
						'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'off'};
					%	m.set.axes = {'xLim', [0, 5], 'xScale', 'log', ...
					%	'yLim', [.04, 10], 'yScale', 'log', 'clipping', 'off'};
				case 'weight' % fit linear data with increasing weight for low ecc.
					m.tasks = 'ratCort fitlm show'; % show regression model
					m.tasks = 'ratCort plot set fitlm pred add'; % fit regression model
					m.x = 'ecc'; m.y = 'fac';
					m.plot.arg = {'o'};
					m.set.axes = {'yLim', [.05, 20], ...
						'xScale', 'log', 'yScale', 'log', 'clipping', 'off'};
					m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};
					m.fitlm.weight = 'weight';
			end
		case 'rat.ecc cum' % ratio of cone count/dens. to ganglion cell count/dens.
			m.tasks = 'densGang countGang densGangSub ratCount plot set';
			m.densGangSub.z = 'count';
			m.ratCount.group = 'type';
			m.x = 'ecc'; m.plot.line = 'type';
			switch 'count'
				case 'count'
					m.y = 'ratio'; m.plot.arg = {'-o'};
					m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
						'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [.3, 20], ...
						'yTick', [.5, 1, 10, 20]};
				case 'dens'
					m.y = 'densDeg'; m.plot.arg = {'o-'};
					m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
						'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4]};
			end
		case 'rat.ecc mid' % midget / all ganglion cells: Grünert (93), Dacey (94)
			m.tasks = 'ratMidget stop plot set';
			m.plot.arg = {'-o'};
			switch 'deg' % unit for eccentricity
				case 'deg' % degrees
					m.x = 'eccDeg'; m.y = 'ratio';
					m.set.axes = {'xLim', [0, 75], 'yScale', 'log', 'yLim', [.1, 1]};
				case 'mm' % (mostly) as published
					m.x = 'ecc'; m.y = 'ratioPer';
					m.set.axes = {'xLim', [0, 15], 'yScale', 'log', 'yLim', [1, 100]};
			end
		case 'resp.freq' % spatial frequency resp.: Croner (95), Wool (18), Yeh (95)
			m.tasks = 'resp plot set'; % plot empirical data
			m.resp.source = 'Croner'; % data source: Croner, Wool, or Yeh
			m.y = 'amp'; % polar component: amp or phase
			m.plot.group = 'source'; m.plot.arg = {'o'};
			switch m.resp.source
				case 'Croner' % Croner (95)
					m.x = 'freqS';
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
						'yScale', 'log' 'yLim', [1, 100]};
				case 'Wool' % Wool (18)
					m.x = 'freqS';
					switch m.y
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yScale', 'log' 'yLim', [1, 100]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yLim', [-180, -90], 'yTick', [-180, -135, -90]};
					end
				case 'Yeh' % Yeh (95)
					m.x = 'freqT';
					switch m.y
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
								'yScale', 'log', 'yLim', [.1, 10]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
								'yLim', [-720, 90], 'yTick', [-540, -360, -180, 0]};
					end
			end
		case 'sens' % ganglion cell contrast sensitivity: Croner (95)
			%	m.tasks = 'respGang'; Wool (18)
			% m.p.kSens ~ (95.788 / .45) / m.p.kRect; cont. sens. (mV/contrast-unit)
			s = 100 * .963; % cont. sens. (Hz/contrast-unit): Croner (95) Figure 9b
			s = s / .75; % convert from descending arm to peak of spatial freq. resp.
			s = s / m.p.kRect; % 17.8 mV/contrast-unit
			fprintf(['Ganglion cell contrast sensitivity, ', ...
				'm.p.kSens = %g mV/contrast-unit\n'], s); % display
			return
		case 'sens beta' % contrast sens. of layer 4CBeta exc. cells: Carandini (97)
			r = [3921009, 52; 3921008, 45; 3821019, 9; 3821021, 28];
				% cell number, response (Hz) at contrast = .3
			m.p.kSensBeta = mean(r(:, 2)) / .3; %  contrast sensitivity
			return
		case 'sep.ecc' % calc. separation of 4CBeta cells from g.c. dens.: pragmatic
			m.tasks = 'sepBeta'; % print results and enter in readTab
			m.x = 'ecc'; m.y = 'sep';
			%	m.sepBeta.mult = 2.5; % 4CBeta cell density is multiple of g.c. dens.
			m.sepBeta.num = 1000; % number of 4CBeta excitatory cells
		case 'tau' % time constant: Ringach (98)
			ss = 6; % number of stages up to, and including, 4CBeta excitatory cells
			tau = .0625; % time to peak of 4CBeta cell (s): Figure 2
			m.p.tau = tau / (ss - 1); % time to peak in a cascade of lowpass filters
			fprintf('Time constant (s): %g\n', m.p.tau); % display
			%	*** This displays .0125. Why is the value .0123 in setLit.m? ***
			return
		case 'wave.x.y' % show preferred wavelengths
			m.tasks = 'select plot set';
			m.x = 'x'; m.y = 'y'; m.z = 'z';
			m.folder = [userpath, '/Data/Prim']; m.p.file = 'Image 5';
			m.select.file = 100;
			m.plot.group = 'name';
			m.plot.funFun = @(d, m) image(flipud(d.wave.wave));
			m.set.axes = {'xLim', [0, 500], 'yLim', [0, 500], 'clipping', 'off'};
			m.set.colorbar = {};
			m.set.title = {'Interpreter', 'none'};
		case 'width.ecc' % determine visual field width for fixed number of cones
			m.tasks = 'densCone select wid plot set'; % *** simpler to use readTab? **
			m.x = 'eccDeg'; m.y = 'width';
			m.select.quad = 'temporal';
			m.wid.cones = 200; % required number of cones
			m.plot.arg = {'-'};
			m.set.axes = {'xLim', 55 * [0, 1], 'yLim', 1.4 * [0, 1]};
	end

	% Set up plotting for publication
	if contains(m.tasks, 'pub') % task list includes pub
		m.pub.units = 'centimeters'; % define unit before using it
		m.pub.axisLength = [5, 5]; % length of x-, y-axes (cm)
		%	m.pub.export = {'colorSpace', 'cmyk'}; % for file export
		m.pub.lineWidth = 1.5; % width of axes and lines (cm)
	end
	
	% Run tasks
	if isfield(m.p, 'file') % m.p.file is defined
		load([m.folder, filesep, m.p.file], 'd'); % load data
	else
		d = table; % empty data table
	end
	stream(d, m); % run tasks

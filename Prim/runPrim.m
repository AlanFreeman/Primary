function runPrim % analyse the behaviour of a visual signal processing model

% Analyse the behaviour of the visual signal processing model, version 'Prim'
% (short for primary visual cortex). This function sets the metadata and
% initiates the stream of analysis tasks.

	% Initialise metadata and read synaptic modulation factors
	import prim.setFig; % find function
	m = runPrimFun; % set function handles
	m = prim.setLit(m); % set metadata calculated from the literature
	m = prim.setDef(m); % set default model parameters
	m = setFig(m, init = []); % set default plotting parameters
	[d, m] = prim.setMod(68, [], m); % read modulation file or set eccentricity

	% Specify the analysis tasks to be performed
	switch 'pc.pc'
		case 'cdi' % plot centre-of-mass displacement index *** to be continued ***
			m.tasks = 'select resp fund crop unpack max index';
			%	'scat desc plot set mark pub'];
			m.x = 'soi'; m.y = 'cdi'; % index: coi, soi, ...
			%	col = ''; % add colour showing index: '' or index name (e.g. eri, osi)
			%	m.mark.loc = [m.p.c.b; m.p.c.c]; % locations to mark: [] or locations
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.cont = m.p.cont;
			m.resp.dir = m.p.c.dir;
			m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = .025;
			m.unpack.var = 'cont';
			m.max.group = 'cont';
			m.index.name = {'ari', 'cdi'}; % doi requires coi and soi
		case 'cor.t' % plot cross-correlation between input and output time courses
			m.tasks = 'select resp corCourse reduce desc plot set pub';
			m.x = 'time'; m.y = 'resp';
			loc = m.p.c.a;
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.p.dir = m.p.c.dirPref; m.p.freqS = m.p.c.freqSPref;
			m.resp.stage = {'genOn', 'betaEx'};
			m.reduce.locO = loc; m.reduce.locI = loc;
			m.plot.funFun = @(d, m)plot(m.p.t - .25, fftshift(d.resp(1, :)), ...
				'clipping', 'on');
			m.set.axes = setFig(m, xSym = .2, ySym = 2);
		case 'count.dom' % plot histogram of cone dominance in centre, surround
			m.tasks = 'resp fund comp crop unpack hist prep plot set';
			m.x = 'resp'; m.y = 'count';
			switch 'cen' % mechanism: cen or sur
				case 'cen', m.resp.stage = 'gangOn'; % centre
				case 'sur', m.resp.stage = 'back'; % surround
			end
			norm = 0; % normalise: 0 or 1
			%	d = table; % empty modulation factor file, for subcortex without cortex
			m.p.wid = .75;
			%	m.p.ecc = 3; m.p.wid = 1;
			m.p.freqS = 0; % uniform stimulus field
			m.p.back = 0; % open loop
			m.resp.cont = .3 * eye(3); % L, M, and S cone contrasts
			m.fund.prop = 'complex'; % amp, complex, phase
			m.crop.z = 'resp';
			m.unpack.var = 'cont';
			m.hist.group = 'cont'; m.hist.z = 'count';
			m.hist.edges = linspace(0, 1, 11);
			m.prep.group = 'cont'; m.prep.z = 'count';
			m.plot.group = 'cont';
			m.plot.funFun = @(d, m)plot(d.(m.x), d.(m.y), ...
				'color', d.cont / max(d.cont));
			x = {'xLim', [0, 1], 'xTick', [0, .5, 1]};
			switch norm % normalise histogram: 0 or 1
				case 0 % no
					m.prep.prep = '';
					m.set.axes = [x, {'yLim', [0, inf]}];
				case 1 % yes
					m.prep.prep = 'norm';
					m.set.axes = [x, {'yLim', [0, 1], 'yTick', [0, .5, 1]}];
			end
		case 'count.ind' % plot histogram of response index
			m.tasks = ...
				'select resp fund crop unpack max index dep hist prep desc plot stop set';
			m.x = 'resp'; m.y = 'count';
			group = ''; % grouping variable: '', 'cycle' *** check cycle option ***
			index = 'maxA'; % coi, soi, ...
			norm = 0; % normalise: 0 or 1
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.cont = m.p.c.cont; % index requires all contrasts
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = .025;
			m.unpack.var = 'cont';
			m.max.group = 'cont';
			m.index.name = index;
			m.dep.name = [index, '.val'];
			m.hist.z = 'count';
			m.set.axes = setFig(m, xZero = 1);
			switch group
				case 'cycle'
					m.select = rmfield(m.select, group);
					m.resp.group = 'cycle';
			end
			if ~ isempty(group)
				g = group;
				m.crop.group = g; m.unpack.group = g; m.max.group = ['cont', g];
				m.index.group = g; m.dep.group = g; m.hist.group = g; m.prep.group = g;
				m.plot.line = g;
			end
			switch norm, case 1 % *** check ***
				m.prep.prep = 'norm'; m.plot.group = 'total';
				m.set.axes = setFig(m, 'yZero', 1);
			end
		case 'count.freqS' % plot histogram of preferred spatial frequency
			m.tasks = ...
				'select resp fund maxOld interp crop hist prep desc plot set';
			m.x = 'freqS';
			switch 'raw' % alternative tasks
				case 'interp' % interpolate on tuning curve
					m.doMax.out = 'freqSTun';					
				case 'raw' % no interpolation
					m.tasks = 'select resp fund doMax crop hist prep desc plot set';
					m.doMax.out = 'freqSPref';
			end
			m.select.phase = d.phase(end); m.select.cycle = m.p.cycle(end);
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.fund.prop = 'amp'; % amp, complex, phase
			m.interp.z = 'freqSTun'; m.interp.out = 'freqSPref';
			m.crop.z = 'freqSPref'; m.crop.radius = 0;
			m.hist.x = 'freqSPref'; m.hist.z = 'count'; m.hist.edges = 10;
			m.prep.prep = 'norm'; m.prep.z = 'count';
			m.plot.group = {'ecc', 'total'};
			m.plot.x = 'freqSPref'; m.plot.y = 'count'; m.plot.z = 'count';
			m.set.axes = {'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
			m.set.axes = {'xLim', [0, 10], 'xTick', [0, 5, 10], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
		case 'count.wid' % plot histogram of direction bandwidth
			m.tasks = 'select resp fund maxOld interp hist desc plot set';
			%	add crop?
			m.x = 'dir';
			m.select.cycle = m.p.cycle(end);
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.doMax.out = 'dirTun';
			m.interp.z = 'dirTun'; m.interp.out = 'band';
			m.hist.x = 'band'; m.hist.z = 'count';
			m.plot.group = {'stage'};
			m.plot.x = 'band'; m.plot.y = 'count'; m.count.z = 'count';
			lim = 70; m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1]};
		case 'dens.x.y' % plot map of cone density
			m.tasks = 'select dens desc plot set';
			m.x = 'x'; m.y = 'y';
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, ...
				shiftdim(d.dens(:, :, :, 2), 1)');
			m.set.axes = setFig(m, xSym = .2, ySym = .2);
		case 'dir.x.y' % plot map of preferred orientation
			m.tasks = 'select resp fund crop max interp prep desc plot set';
			ana = ''; % analysis variant: '', 'dev'
			map = 'rgb'; % mapping from direction to colour code: hsv, rgb
			m.select.cycle = m.p.cycle(end);
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = .005;
			m.max.out = 'dirTun';
			m.interp.x = 'dir'; m.interp.z = 'dirTun'; m.interp.out = 'dirPref';
			m.prep.x = 'loc'; m.prep.z = 'dirPref';
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.dirPref));
			w = .5 * m.p.wid; lim = w * [-1, 1]; tick = w * [-1, 0, 1];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colormap = 'hsv';
			switch ana % analysis variant
				case 'dev' % development
					m.select = rmfield(m.select, 'cycle');
					g = 'cycle';
					m.resp.group = g; m.crop.group = g; m.doMax.group = g;
					m.interp.group = g; m.prep.group = g; m.plot.group = g;
			end
			switch map % mapping from direction to colour code
				case 'hsv' % response strength coded by saturation
					m.max.struct = 1; % doMax has multiple outputs
					m.prep.prep = {'sat', 'image'}; m.prep.var = 'dir';
					m.set.colorbar = {'xLim', [0, 1], 'ticks', [0, .5, 1], ...
						'tickLabels', {'-90', '0', '90'}};
				case 'rgb' % response strength uncoded
					m.prep.prep = {'orient', 'image'};
					m.set.caxis = 90 * [-1, 1];
					m.set.colorbar = {'xLim', 90 * [-1, 1], 'ticks', 90 * [-1, 0, 1]};
			end
		case 'dir.x.y pred' % map of pref. orient. predicted from g.c coverage * fix
			m.tasks = 'select dens predOr desc plot set';
			m.select.cycle = m.p.cycle(end);
			%	m.dens.rad = .07;
			m.orient.z = 'dens';
			%	m.orient.debug = 1; m.p.dir = -36; m.p.freqS = 4.9;
			dir = -90: 1: 90; m.orient.dir = dir(1: end - 1);
			m.orient.freqS = 12;
			m.x = 'x'; m.y = 'y'; m.z = 'orPref'; % m.z = 'freqSPref';
			m.plot.group = {'ecc', 'freqS'};
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.(m.z), 1));
			w = .25;
			lim = w * [-1, 1]; tick = w * [-1, 0, 1];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			%	m.set.caxis = [-1.5, 1.5]; m.set.colorbar = {'xLim', [-1.5, 1.5]};
			m.set.colorbar = {'xLim', 90 * [-1, 1], 'ticks', 90 * [-1, 0, 1]};
			m.set.colormap = 'hsv'; % m.set.caxis = limC * [-1, 1];
		case 'fig' % plot figureind
			switch 'conv'
				case 'conv' % convergence functions
					m.tasks = 'fig plot set pub';
					m.x = 'x'; m.y = 'f';
					m.fig.conv = 1; m.fig.xLim = .1;
					 m.plot.group = 'source';
					m.set.axes = setFig(m, xSym = .2, yZero = 1);
				case 'course' % stimulus time course *** fix ***
					m.tasks = 'stim fig plot set';
					m.fig.course = 1;
					m.x = 't'; m.y = 'course';
					yLim = m.p.cont(1);
					m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
						'yLim', yLim * [-1, 1], 'yTick', yLim * [-1, 0, 1]};
			end
		case 'ind.ind' % plot response index versus response index
			m.tasks = ['select resp fund crop unpack max index ', ...
				'scat desc plot set pub'];
			m.x = 'soi'; m.y = 'coi'; % index: coi ...
			col = ''; % add colour showing index: '' or index name (e.g. emi)
			mark = 1; % mark locations: '' or 1
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.cont = m.p.c.cont; m.resp.dir = m.p.c.dir;
			m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx'; m.resp.stage = 'gangOff';
			m.crop.radius = .025;
			m.unpack.var = 'cont';
			m.max.group = 'cont';
			m.index.name = {m.x, m.y, col}; % doi requires coi and soi
			m.scat.name = m.index.name;
			m.plot.funFun = @(d, m) scatter(d.(m.x), d.(m.y), [], d.(col), ...
				'filled', markerEdgeColor = 'k', lineWidth = .25, userData = d.loc);
			m.set.axes = setFig(m, xZero = 1, yZero = 1);
			m.set.colorbar = {'title', col};
			%	z = setFig(m, figIndexB = 8);
			%	m.set.caxis = z{1}; m.set.colorbar = [m.set.colorbar, z{2: 3}];
			if isempty(col) % no colour showing index value
				m.index.name = {m.x, m.y};
				switch col % which index?
					case 'cdi' % *** fix ***
						m.plot.funFun = @(d, m) scatter(d.resp.(m.x).val, ...
							d.resp.(m.y).val, [], d.resp.cdi.disp.val, 'filled');
					otherwise
						m.plot.funFun = @(d, m) scatter(d.(m.x), d.(m.y), ...
							'filled', markerEdgeColor = 'k', userData = d.loc);
				end
				m.set = rmfield(m.set, 'colorbar');
			end
			if ~ isempty(mark) % mark locations
				m.tasks = replace(m.tasks, 'set', 'set mark');
				m.mark.loc = m.p.c.e; % locations to mark: ls x 2
				%	m.mark.loc = [m.p.c.b; -.048, .048; .024, -.016]; % double opponent
			end
		case 'ind.x.y' % plot map of response index
			m.tasks = ['select resp fund crop unpack max index pca ', ...
				'dep prep desc plot set'];
			m.x = 'loc'; m.y = 'loc';
			col = 'pc1'; % add colour to each point: index name
			name = {'eti', 'mti'}; % names of indices to include in PCA
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.cont = m.p.c.cont;
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = .025;
			m.unpack.var = 'cont';
			m.max.group = 'cont';
			m.index.name = [name, {'pca', col}];
			m.pca.name = name; m.pca.reverse = 1; m.pca.method = 'pca';
			m.dep.name = [col, '.val'];
			m.prep.prep = 'image';
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.resp, 1));
			lim = .2; m.set.axes = setFig(m, xSym = lim, ySym = lim);
			%	m.set.caxis = lim;
			m.set.colorbar = {'title', col};
			switch col % cdi is a special case
				case 'cdi', m.dep.name = 'cdi.val.disp';
			end
		case 'loc' % find the location of a cell marked by a datatip
			m.tasks = 'loc';
		case 'max.x.y' % plot map of maximal response or selectivity index
			m.tasks = 'select resp fund crop unpack max prep desc plot set';
			m.x = 'loc'; m.y = 'loc'; m.z = 'resp';
			ana = ''; % analysis type: '', chrom, cont, image, orient *** fix ***
			desat = 0; % 1 to desaturate poor orientation selectivity, otherwise 0
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = 0; m.crop.radius = .02;
			m.prep.prep = 'image';
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.resp));
			m.set.axes = setFig(m, xSym = m.p.wid, ySym = m.p.wid);
			m.set.colorbar = {}; %	m.set.colorbar = {'ticks', lim * [0, .5, 1]};
			switch ana
				case 'chrom' % calculate chromaticity selectivity index
					m.tasks = replace(m.tasks, 'max', 'max chrom'); % add chrom task
					m.resp.cont = m.p.contMag * [1, -1, 0; 1, 1, 1];
					m.unpack.var = 'cont';
					m.max.group = 'cont';
				case 'cont' % calculate maximum for achromatic and equiluminant stimuli
					m.resp.cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1];
					m.unpack.var = 'cont';
					m.max.group = 'cont'; m.prep.group = 'cont';
					m.plot.group = {'ecc', 'cont'};
				case 'image' % *** doesn't work because max doesn't recognise files ***
					m.p.stimS = 'image'; m.p.time = .25;
					m.resp.file = 'Image 8'; % image library
					m.resp.index = (2)'; % image file numbers
					m.resp = rmfield(m.resp, 'freqS'); % not relevant
					m.resp.win = .75; % image window (fraction of image width)
				case 'orient' % calculate orientation selectivity index
					m.tasks = replace(m.tasks, 'max', 'max orient'); % add orient task
					m.max.out = 'dirTun';
			end
			switch desat % code OSI by saturation in plot
				case 1
					m.max.struct = 1; % max outputs both maximum and tuning curve
					m.prep.prep = {'sat', 'image'}; m.prep.var = 'max';
			end
		case 'max.x.y image' % plot map of maximal response to images
			m.tasks = 'select resp peak crop max dep prep desc plot set';
			m.x = 'loc'; m.y = 'loc'; m.z = 'resp';
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.p.stimS = 'image'; m.p.time = .25;
			m.resp.file = 'Image 8'; % image library
			m.resp.index = (1: 5)'; % image index numbers
			m.resp.index = [1; 2];
			m.resp.dir = m.p.c.dir; % image drift direction (deg)
			m.resp.dir = [144; 162];
			%	m.resp.win = .75; % image window (fraction of image width)
			m.resp.domain = 'time'; m.resp.stage = 'betaEx';
			m.crop.radius = 0; m.crop.radius = .02;
			%	m.dep.name = 'max';
			m.dep.name = 'max.val';
			m.prep.prep = 'image';
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.resp));
			m.set.axes = setFig(m, xSym = m.p.wid, ySym = m.p.wid);
			m.set.colorbar = {}; % {'ticks', [-180, 0, 180]};
		case 'mod' % calculate and save synaptic modulations
			m.tasks = 'mod save';
			m.mod.serial = 74; % serial number of output file
			m.save.name = ['Modulation ', num2str(m.mod.serial), '.mat'];
			switch 2 % development period
				case 1 % drifting achromatic gratings, simulating previsual dev't
					m.mod.dir = m.p.c.dir; m.mod.freqS = m.p.c.freqS;
					m.mod.init = 'init'; % initialise output file
				case 2 % drifting equiluminant gratings, simulating quasivisual dev't
					m.mod.cont = m.p.contMag * [1, -1, 0];
					m.mod.dir = m.p.c.dir; m.mod.freqS = m.p.c.freqS;
					m.mod.init = 'cont'; % cont or init
				case 3 % pulsed natural images, simulating visual development
					m.p.stimS = 'image'; m.p.stimT = 'pulse'; m.p.time = .25;
					m.mod.init = 'cont'; % cont or init
					m.mod.file = 'Image 8';
					m.mod.files = []; % index numbers of image files
					%	m.mod.save = 1; % save single-cycle files
			end
		case 'mod.x.y dom' % plot synaptic mod'n map, colour with cone dominance
			m.tasks = 'resp fund comp save'; % save cone dominance
			switch 'ungroup'
				case 'off'
					m.resp.stage = 'gangOff';
					m.save.name = 'ConeDomOff';
				case 'on'
					m.resp.stage = 'gangOn';
					m.save.name = 'ConeDomOn';
				case 'plot'
					m.tasks = 'select reduce unpack colMod plot set add';
					m.p.file = 'Modulation.mat';
					m.plot.group = 'stage';
				case 'ungroup'
					m.tasks = 'select reduce unpack colMod plot set add';
					m.p.file = 'Modulation.mat';
			end
			m.p.ecc = 3;
			m.p.freqS = 0; m.p.back = 0; m.resp.cont = .3 * eye(3);
			m.select.cycle = 3;
			m.fund.prop = 'complex'; % amp, complex, phase
			m.reduce.locBeta = [0, 0]; m.reduce.locBeta = [-.093, -.074];
			m.unpack.var = 'locGen'; m.unpack.locGen = {'x', 'y'};
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @ (d, m) scatter(d.x, d.y, max(1, 36 * d.resp), ...
				d.colour, 'filled');
			w = .5 * m.p.wid; lim = w * [-1, 1]; tick = w * [-1, 0, 1];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.add.funFun = @(d, m) plot(m.reduce.locBeta(1), m.reduce.locBeta(2), ...
				'ok', 'markerFaceColor', 'k', 'markerSize', 10);
		case 'op.x.y' % plot map of opponency index
			m.tasks = 'load reduce unpack field op prep desc plot set';
			m.load.name = 'Resp pulse 56.mat';
			m.reduce.cont = [0, .3, 0; .3, 0, 0];
			m.unpack.var = 'cont'; m.field.group = 'cont';
			m.prep.prep = 'image'; m.prep.x = 'loc'; m.prep.z = 'oi';
			m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.oi));
			m.set.axes = setFig(m, xSym = m.p.wid, ySym = m.p.wid);
			m.set.colorbar = {};
		case 'pc.pc' % principal components analysis
			m.tasks = ['select resp fund crop unpack max index pca ' ...
				'scat desc plot set'];
			m.x = 'pc1'; m.y = 'pc2'; %	m.y = 'pc3';
			col = 'soi'; % add colour showing index: index name
			mark = 0; % mark locations: '' (default is no marked locations) or 1
			method = ''; % analysis method: '' (default is pca) or 'tsne'
			name = {'eti', 'mti'}; % names of indices to include in PCA: original
			%	name = {'eti', 'lti', 'mti'}; % arc over and under
			%	name = {'ati', 'eti', 'lti', 'mti'}; % right angle, big pc3
			%	name = {'ami', 'emi', 'ati', 'eti', 'lti', 'mti'};
				% names of indices to include in PCA
			type = ''; % plot type: '' (default is pc2 v. pc1), biplot, scree
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.cont = m.p.c.cont; m.resp.dir = m.p.c.dir;
			m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.crop.radius = .025;
			m.unpack.var = 'cont';
			m.max.group = 'cont';
			m.index.name = [name, {'pca', col}];
			m.pca.name = name; m.pca.reverse = 1;
			m.scat.name = {m.x, m.y, col}; % m.scat.descend = 1;
			m.plot.funFun = @(d, m) scatter(d.(m.x), d.(m.y), [], d.(col), ...
				'filled', markerEdgeColor = 'k', userData = d.loc, clipping = 'off');
			m.set.axes = setFig(m, xZero = 1, yZero = 1);
			m.set.colorbar = {'title', col};
			switch col % set colour limits
				case 'osiA', m.set.caxis = [0, .6];
			end
			switch method % analysis method
				case 'tsne'
					m.pca.method = 'tsne';
					m.set = rmfield(m.set, 'axes');
				otherwise, m.pca.method = 'pca';
			end
			switch mark % mark a location on the map
				case 1
					m.tasks = replace(m.tasks, 'set', 'set mark');
					m.mark.loc = [m.p.c.b; m.p.c.c; m.p.c.d]; % locations to mark
			end
			switch type % plot type: '', biplot, scree
				case 'biplot'
					label = ["0", "2", "4", "6", "8", "10", "12", "14", "16", "18"];
					m.plot.funFun = @(d, m) biplot(d.resp.vec.val(:, 1: 2), ...
						varLabels = label);
					m.set.axes = {'xLim', .5 * [-1, 1], 'yLim', [-.1, .5]};
				case 'scree'
					m.tasks = replace(m.tasks, 'scat ', '');
					m.plot.funFun = @(d, m) pareto(d.resp.frac.val);
					m.set.axes = setFig(m, yZero = 100);
					m.set = rmfield(m.set, 'colorbar');
			end
		case 'rad.ecc' % calculate radius from spatial frequency response
			m.tasks = 'resp reduce fund radius save list'; % save radius data
			m.tasks = 'load plot set'; % plot it
			m.x = 'ecc'; m.y = 'radius'; m.z = 'radius';
			m.p.file = '';
			m.p.ecc = 30; % save one eccentricity at a time
			m.p.back = 0;
			switch m.p.ecc
				case .1, m.p.wid = .2; xLim = log10([2, 50]);
				case .3, m.p.wid = .2; xLim = log10([2, 50]);
				case 1, m.p.wid = .5; xLim = log10([2, 50]);
				case 3, m.p.wid = .5; xLim = log10([2, 50]);
				case 10, m.p.wid = 1; xLim = log10([.5, 12]);
				case 30, m.p.wid = 2.5; xLim = log10([.1, 10]);
			end
			switch 'cen' % choose mechanism
				case 'cen' % centre mechanism
					m.resp.stage = 'gangOff'; % centre
					yLim = [.007, .3]; yTick = [.01, .1];
				case 'sur' % surround mechanism
					m.resp.stage = 'back'; % surround
					yLim = [.007, 5]; yTick = [.01, .1, 1];
			end
			freqS = logspace(xLim(1), xLim(2), 30); freqS = [0, freqS];
			m.resp.freqS = freqS';
			m.reduce.loc = [0, 0];
			m.fund.prop = 'amp'; % amp, complex, phase		
			m.save.name = 'Radius.mat'; m.load.name = m.save.name;
			m.plot.arg = {'o'};
			m.plot.group = 'stage';
			m.set.axes = {'xScale', 'log', 'xLim', [.006, 40], ...
				'xTick', [.01, .1, 1, 10], 'yScale', 'log', 'yLim', yLim, ...
				'yTick', yTick, 'clipping', 'off'};
		case 'ratio.x.y' % plot map of L-cone component in centre and surround
			m.tasks = 'resp fund comp reduce unpack plot set';
			m.tasks = 'resp fund comp reduce crop unpack plot set';
			switch 'cen' % mechanism
				case 'cen' % centre
					m.p.arrayC = 'gangOn'; m.reduce.stage = 'gangOn';
				case 'sur' % surround
					m.p.arrayC = 'cone'; m.reduce.stage = 'back';
			end
			m.p.freqS = 0; % uniform stimulus field
			m.p.kSur = 0; % open loop
			m.resp.cont = .3 * eye(3); % all cone contrasts
			m.fund.prop = 'complex'; % amp, complex, phase
			m.reduce.cont = .3 * [1, 0, 0]; % L-cone only
			m.crop.z = 'resp';
			m.unpack.var = 'loc'; m.unpack.loc = {'x', 'y'};
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)scatter(d.x, d.y, 60 * d.resp, ...
				d.resp .* [1, 0, 0], 'filled');
			w = .5 * m.p.wid; % half-width of visual field (deg)
			lim = w * [-1, 1]; tick = w * [-1, 0, 1];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
		case 'resp' % save file for excitatory cell response to pulsed stimulus
			m.tasks = 'select resp peak save'; 
			m.p.stimS = 'grating'; % m.p.stimS = 'image';
			m.p.stimT = 'pulse'; m.p.time = .25;
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.domain = 'time'; m.resp.stage = 'betaEx';
			switch m.p.stimS
				case 'grating'
					m.resp.cont = m.p.c.cont;
					dir = linspace(-90, 90, 20 + 1); m.resp.dir = dir(1: end - 1)';
					freq = m.p.c.freqS; m.resp.freqS = freq(2: end); % remove 0 cyc/deg
					m.resp.phaseS = [-180, -90, 0, 90]';
					m.save.name = 'Resp pulse 68.mat';
				case 'image'
					m.resp.file = 'Image 2';
					m.resp.image = 1: 1534; % m.resp.image = 1: 100;
					m.save.name = 'Field image solveC.mat';
			end
		case 'resp.dir' % plot direction tuning
			m.tasks = 'select resp reduce fund unpack max desc plot set pub';
			m.x = 'dir'; m.y = 'resp'; m.z = 'resp';
			ana = 'dev'; % analysis type: '', cont, dev, devEqui, fit
			loc = m.p.c.b;
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.reduce.loc = loc;
			m.plot.group = {'loc', 'stage'};
			m.plot.funFun = @(d, m)plot(d.dir, d.resp.dir.tun.val);
			yLim = 4; m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', yLim * [0, 1], 'yTick', 4 * [0, .5, 1], 'clipping', 'off'};
			switch ana
				case 'cont'
					m.resp.cont = m.p.c.cont;
					m.unpack.var = 'cont';
					m.max.group = 'cont';
					m.plot.group = {'loc', 'stage'}; m.plot.line = 'cont';
				case {'dev', 'devEqui'} % development over cycles
					m.select = rmfield(m.select, 'cycle');
					switch ana, case 'devEqui', m.p.cont = m.p.contMag * [1, -1, 0]; end
					m.resp.group = 'cycle'; m.reduce.group = 'cycle';
					m.max.group = 'cycle'; m.plot.line = 'cycle';
				case 'fit' % add interpolated fit
					m.tasks = [m.tasks, ' interp add'];
					m.interp.z = 'dirTun';
					m.add.y = 'dirTun';
					m.add.funFun = @(d, m)plot(d.(m.x), squeeze(d.dirTun));
					m.plot.funFun = @(d, m)plot(d.(m.x), squeeze(d.dirTun), 'o');
			end
		case 'resp.dir.freqS' % plot direction and spatial frequency tuning
			m.tasks = 'select resp reduce fund unpack desc plot set';
			if isempty(m.p.file), m.tasks = erase(m.tasks, 'select '); end
			m.x = 'dir'; m.y = 'freqS'; m.z = 'resp';
			ana = 'cont'; % analysis type: '', cont *** group on phaseS ***
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.dir = linspace(-180, 180, 20 + 1)'; m.resp.freqS = m.p.c.freqS;
			m.resp.stage = 'betaEx';
			m.reduce.loc = m.p.c.b;
			m.unpack.var = 'loc';
			m.plot.group = {'loc', 'stage'};
			m.plot.funFun = @(d, m)imagesc(d.dir, d.freqS, squeeze(d.resp)');
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', m.p.c.freqSRange}; % 'yTick', [0, 5, 10]};
			m.set.colorbar = {};
			switch ana
				case 'cont'
					m.resp.cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1];
					m.unpack.var = {'cont', 'loc'};
					m.plot.group = {'loc', 'stage', 'cont'};
			end
		case 'resp.freqS' % plot spatial frequency tuning with optimal direction
			m.tasks = 'select resp reduce fund unpack max desc plot set pub';
			ana = 'cont'; % analysis type: '', 'cont', 'cycle'
			gang = 1; % 0 for cortex, 1 for ganglion cell
			m.select.phase = d.phase(end); m.select.cycle = m.p.cycle(end);
			m.p.cont = m.p.contMag * [1, -1, 0]; % equiluminance
			m.resp.dir = m.p.c.dir; m.resp.freqS = m.p.c.freqS;
			f = [1, 20]; f = log10(f); m.resp.freqS = logspace(f(1), f(2), 20)';
			m.resp.stage = 'betaEx'; m.resp.stage = 'gangOff';
			m.reduce.loc = m.p.c.e;
			m.plot.group = {'loc', 'stage'};
			m.plot.funFun = @(d, m)plot(d.freqS, d.resp.freqS.tun.val');
			x = {'xScale', 'log', 'xLim', [1, 20], 'xTick', [1, 10]};
			y = {'yScale', 'log', 'yLim', [.1, 7], 'yTick', [.1, 1, 10]}; % Index
			%	y = {'yScale', 'log', 'yLim', [.4, 11], 'yTick', [.1, 1, 10]}; % DevPost
			%	y = {'yScale', 'log', 'yLim', [.5, 11], 'yTick', [.1, 1, 10]}; % PCA
			m.set.axes = [x, y];
			switch ana % alternative tasks
				case 'cont'
					m.resp.cont = m.p.c.cont;
					m.unpack.var = 'cont';
					m.max.group = 'cont';
					m.plot.group = {'loc', 'stage', 'cont'}; % m.plot.line = 'cont';
				case 'cycle'
					g = 'cycle';
					m.select = rmfield(m.select, g); % all cycles
					m.resp.group = g; m.reduce.group = g; m.max.group = g;
					m.plot.line = g;
			end
			switch gang % ganglion cell or cortex?
				case 1 % ganglion cell
					m.resp.stage = 'gangOff';
					m.set.axes = [x, {'yScale', 'log', 'yLim', [13, 30], ...
						'yTick', [15, 20, 25, 30]}]; % equiluminant, achromatic
					m.set.axes = [x, {'yScale', 'log', 'yLim', [.4, 30], ...
						'yTick', [.1, 1, 10, 30]}]; % L-specific, M-specific
			end
		case 'resp.freqS fix' % plot spatial frequency tuning with fixed directions
			m.tasks = 'select resp reduce fund prep unpack desc plot set'; % uncentred
			%	m.tasks = 'select resp reduce centre fund prep unpack desc plot set';
			m.x = 'freqS'; m.y = 'resp'; m.z = 'resp';
			area = ''; % '', 'sub' for subcortex alone
			m.select.phase = d.phase(end); m.select.cycle = m.p.cycle(end);
			%	m.p.act = 0;
			m.resp.cont = m.p.contMag * [1, -1, 0];
			m.resp.dir = (0: 15: 90)'; m.resp.dir = 0;
			m.resp.stage = {'betaEx'};
			m.reduce.loc = [0, 0];
			m.reduce.loc = m.p.c.big;
			m.fund.prop = 'amp'; % amp, complex, phase
			xLim = 30; freqS = linspace(0, xLim, 50); m.resp.freqS = freqS';
			%	m.resp.back = [0; 1]; % m.resp.back = 0; %	m.resp.stage = 'back';
			m.unpack.var = {'dir', 'loc', 'stage'};
			m.plot.group = 'stage';
			m.plot.line = {'dir', 'loc'};
			m.plot.funFun = @(d, m)plot(d.freqS, squeeze(d.resp)');
			a = {'xLim', xLim * [0, 1], 'xTick', xLim * [0, .5, 1]};
			switch area % subcortical or cortex
				case 'sub'
					%	set Modulation file to [], and set eccentricity in call to setMod
					m.tasks = erase(m.tasks, 'select '); % no need to select title
			end
			switch m.fund.prop % set axes
				case 'amp'
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 13], ...
						'yTick', [1, 3, 10]];
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 8], 'yTick', [1, 3, 8]];
				case 'phase'
					m.prep.unwrap = 'freqS';
					m.set.axes = [a, 'yLim', 200 * [-1, 1], ...
						'yTick', 180 * [-1, 0, 1], 'clipping', 'off'];
					m.set.axes = [a, 'yLim', [120, 180], ...
						'yTick', [120, 150, 180], 'clipping', 'off'];
			end
		case 'resp.freqS multi' % plot spatial frequency tuning, multiple arrays
			m.tasks = 'select resp select centre fund prepOld desc plot set';
				% *** fix centre: location dimensions ***
			m.tasks = 'select resp select fund prepOld desc plot set';
			m.select.cycle = m.p.cycle(end);
			m.p.dir = 132;
			xLim = [1, 20]; % subcortex: [1, 13]; cortex: [1, 20]
			x = log10(xLim); % limits for log x
			freqS = logspace(x(1), x(2), 20); m.resp.freqS = freqS';
			%	m.resp.back = [0; 1]; %	m.resp.back = 0;
			%	m.resp.locS = [-.02; .02];
			%	m.resp.stage = 'back'; m.resp.stage = {'gangOff', 'gangOn'};
			m.resp.stage = {'genOff', 'betaEx'};
			m.select.group = 'stage'; m.select.loc = [0, 0];
			%	m.reduce.loc = [.03, 0; .03, .05; -.04, .05; ...
			%	-.06, 0; -.03, -.05; .03, -.04];
			m.centre.group = {'stage', 'loc'};
			m.fund.group = {'stage', 'loc'}; m.fund.prop = 'amp'; % amp, comp'x, phase
			m.x = 'freqS'; m.y = 'resp';
			m.plot.group = 'stage';
			m.plot.line = {'loc'};
			m.plot.funFun = @(d, m)plot(d.freqS, squeeze(d.resp)');
			switch 'cort' % subcortex or cortex?
				case 'cort', a = {'xScale', 'log', 'xLim', [1, 20], 'xTick', [1, 20]};
				case 'sub', a = {'xScale', 'log', 'xLim', [1, 13], 'xTick', [1, 3, 10]};
			end
			switch m.fund.prop % set axes
				case 'amp'
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 13], ...
						'yTick', [1, 3, 10]];
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 6], 'yTick', [1, 3, 6]];
				case 'phase'
					m.prep.unwrap = 'freqS';
					m.set.axes = [a, 'yLim', 200 * [-1, 1], ...
						'yTick', 180 * [-1, 0, 1], 'clipping', 'off'];
					m.set.axes = [a, 'yLim', [90, 200], ...
						'yTick', [90, 135, 180], 'clipping', 'off'];
			end
		case 'resp.freqS fit' % fit model's spatial frequency response with Gaussian
			m.tasks = 'resp reduce fund prepOld unpack desc plot set';
			%	m.tasks = 'resp reduce fund prepOld unpack desc fitGauss show';
			m.tasks = 'resp reduce fund prepOld unpack desc plot set fitGauss pred add';
			m.p.arrayC = 'gangOn'; m.reduce.stage = 'gangOn';
			switch 'cen' % response component: center, surround, both or empty
				case 'both' % centre and surround components
					m.fitGauss.comp = 'both'; m.fitGauss.coef = 4;
				case 'cen' % centre component
					m.p.back = 0;
					m.fitGauss.comp = 'cen'; m.fitGauss.coef = 1;
				case 'sur' % surround component *** check: task back removed ***
					t = 'resp reduce fund prepOld unpack desc';
					m.tasks = [t,  ' plot set'];
					%	m.tasks = [t, ' fitGauss show'];
					%	m.tasks = [t, ' plot set fitGauss pred add'];
					m.p.arrayC = 'cone'; m.reduce.stage = 'back';
					m.p.back = 0;
					m.p.solver = 'solveF'; % only works for frequency domain
					m.fitGauss.comp = 'sur'; m.fitGauss.coef = 1;
			end
			%	m.p.cont = .3 * [1, 0, 0]; m.p.cont = .3 * [0, 1, 0];
			xLim = log10([.1, 50]); % limits for log x 
			freqS = logspace(xLim(1), xLim(2), 30); m.resp.freqS = freqS';
			%	m.resp.back = [0; 1];
			m.reduce.loc = [0, 0]; % m.reduce.loc = [-.1, .1];
			m.fund.prop = 'amp'; % amp, complex, phase
			m.unpack.var = {'loc', 'stage', 'freqS'};
			%	m.unpack.var = {'loc', 'stage', 'freqS', 'back'};
			m.x = 'freqS'; m.y = 'resp';
			m.plot.group = 'stage'; %	m.plot.line = 'stage';
			%	m.plot.line = 'back';
			m.plot.funFun = @(d, m)plot(d.freqS, squeeze(d.resp));
			a = {'xScale', 'log', 'xLim', 10 .^ xLim};
			m.add.arg = {'o'};
			switch m.fund.prop % set axes
				case 'amp'
					m.set.axes = [a, 'yScale', 'log'];
					m.set.axes = [a, 'yScale', 'log', 'yLim', [.5, 10]];
				case 'phase'
					m.prep.unwrap = 'freqS';
					yLim = 200; m.set.axes = [a, 'yLim', yLim * [-1, 1], ...
						'yTick', yLim * [-1, 0, 1], 'clipping', 'off'];
					%	m.set.axes = a;
			end
		case 'resp.freqT' % plot temporal frequency response
			m.tasks = 'select resp reduce fund unpack desc plot set';
			m.select.cycle = m.p.cycle(end);
			m.p.dir = 60; m.p.freqS = 6.6;
			m.resp.stage = 'betaEx';
			xLim = log10([1, 50]); % limits for log frequency
			freqT = logspace(xLim(1), xLim(2), 30); m.resp.freqT = freqT';
			m.reduce.loc = [0, 0];
			m.fund.prop = 'amp'; % amp, complex, phase
			m.unpack.var = 'stage';
			m.x = 'freqT'; m.y = 'resp';
			m.plot.group = 'stage';
			m.plot.funFun = @(d, m)plot(d.freqT, squeeze(d.resp));
			m.set.axes = {'xScale', 'log'}; % 'yLim', [0, 6]};
		case 'resp.t' % plot model time course
			m.tasks = 'select resp reduce unpack prep desc plot set';
			if isempty(m.p.file), m.tasks = erase(m.tasks, 'select '); end
			m.x = 'time'; m.y = 'resp';
			ana = ''; % analysis type: '', cont, dev, image, phaseS
			m.p.stimS = 'grating'; % m.p.stimS = 'image';
			m.p.stimT = 'drift'; m.p.stimT = 'pulse';
			m.p.solver = 'solveF'; m.p.solver = 'solveT';
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.resp.domain = 'time';
			m.resp.stage = 'betaEx'; % m.resp.stage = 'genOn';
			m.reduce.loc = m.p.c.b;
			m.unpack.var = {'stage', 'loc'};
			m.prep.prep = 'real';
			m.plot.group = 'loc'; m.plot.line = 'stage';
			m.plot.funFun = @(d, m)plot(d.time, squeeze(d.resp));
			switch ana
				case 'cont' % *** direction, frequency setting will match only one cont.
					m.resp.cont = m.p.contMag * [0, 1, 0; 1, -1, 0; 1, 0, 0; 1, 1, 1];
					m.unpack.var = {'stage', 'loc', 'cont'};
					m.prep.group = 'cont';
					m.plot.group = {'cont', 'loc'};
				case 'dev' % development over cycles
					m.select.cycle = m.p.cycle(1: end);
					m.resp.group = 'cycle'; m.reduce.group = 'cycle';
					m.unpack.group = 'cycle';
					m.plot.group = {'loc', 'stage'}; m.plot.line = 'cycle';
				case 'image' % multiple stimuli
					var = {'file', 'dir'}; % stimulus variable
					m.unpack.var = ['stage', 'loc', var]; m.prep.group = var;
					m.plot.group = 'file'; m.plot.line = 'dir';
				case 'phaseS'
					m.resp.phaseS = [-180, -90, 0, 90]';
					m.unpack.var = {'loc', 'phaseS'};
			end
			switch m.p.stimS
				case 'grating'
					m.p.cont = m.p.contMag * [1, -1, 0];
					m.p.dir = m.p.c.dirPref; m.p.freqS = m.p.c.freqSPref;
					m.p.phaseS = 0;
					switch m.p.stimT
						case 'drift'
							switch m.p.solver
								case 'solveF', off = 0;
								case 'solveT', m.p.time = 1; off = .5;
							end
							x = {'xLim', off + .5 * [0, 1], 'xTick', off + .5 * [0, .5, 1]};
							lim = 25; tick = lim;
							if ~ m.p.act % potential
								y = {'yLim', lim * [-1, 1], 'yTick', lim * [-1, 0, 1], ...
									'clipping', 'on'};
							else % impulse rate
								y = {'yLim', lim * [0, 1], 'yTick', tick * [0, .5, 1], ...
									'clipping', 'on'};
							end
							m.set.axes = [x, y]; % m.set.axes = x;
					end
				case 'image'
					m.p.time = .25;
					m.resp.file = 'Image 8'; % image library
					m. resp.index = (2)'; % image file numbers
					m.resp.dir = (0)'; % image drift direction (deg)
					m.resp.win = .75; % image window (fraction of image width)
					m.set.axes = setFig(m, xZero = m.p.time);
			end
		case 'resp.resp' % polar plot g.c. response to S-stim and lum-stim
			m.tasks = 'resp fund crop centre unpack plot set export';
			m.p.file = '';
			%	m.p.ecc = 3; m.p.wid = .35; rLim = [0, 23]; rTick = 10;
			%	m.p.ecc = 10; m.p.wid = .5;
			m.p.ecc = 30; m.p.wid = 2; rLim =[0, 12]; rTick = 5;
			m.p.freqS = 0;
			m.resp.stage = 'gangOff';
			m.resp.cont = m.p.contMag * [0, 0, 1; 1, 1, 0]; % S-, luminance-stimulation
			m.fund.prop = 'complex'; % amp, complex, phase
			m.crop.z = 'resp';
			m.unpack.var = 'cont';
			m.x = 'resp'; m.y = 'resp'; m.z = 'resp';
			m.plot.line = {'stage', 'cont'};
			m.plot.funFun = @(d, m)polarplot(squeeze(d.resp), 'o');
			m.set.axes = {'rLim', rLim, 'rTick', rTick * [1, 2]};
		case 'resp.resp cart' % g.c. response to S-stim vs lum-stim response
			m.tasks = 'resp fund crop unpack plot set';
			m.p.ecc = 3; m.p.wid = .35;
			m.p.ecc = 10; m.p.wid = .5;
			m.p.freqS = 0;
			m.resp.stage = 'gangOff';
			m.resp.cont = .3 * [0, 0, 1; 1, 1, 0]; % S-, luminance-stimulation
			m.fund.prop = 'amp';
			m.crop.z = 'resp';
			m.unpack.var = 'cont';
			m.x = 'resp'; m.y = 'resp'; m.plot.group = 'stage';
			m.plot.funFun = @(d, m)plot(squeeze(d.resp(2, 1, :)), ...
				squeeze(d.resp(1, 1, :)), 'o');
			m.set.axes = {'xLim', [0, 6], 'xTick', [0, 3, 6], ...
				'yLim', [0, 6], 'yTick', [0, 3, 6]};
		case 'resp.x.y' % plot receptive field: cross corr'n of stim., response
			m.tasks = 'load reduce unpack field prep desc plot set add';
			loc = m.p.c.b; % cell location (deg)
			span = 1; % 1 for standard view, .6 for expanded view around cell
			type = 'mono'; % plot type: 'contour', 'mono', 'true'
			m.load.name = 'Resp pulse 67.mat';
			m.reduce.loc = loc;
			m.unpack.var = 'cont'; m.field.group = 'cont';
			m.prep.group = 'cont'; m.prep.prep = 'image'; m.prep.x = 'locS';
			m.plot.group = 'cont';
			m.set.axes = setFig(m, xSym = .2 + 1e-6, ySym = .2);
				% xSym = .2 gives PDF with uniform image; why ??
			m.add.funFun = @(d, m) plot(loc(1), loc(2), 'ok', LineWidth = 2);
			%	m.pub.export = {};
			if span < 1 % expanded view around loc
					span = span * .5 * m.p.wid * [-1, 1];
					lim = {'xLim', loc(1) + span, 'yLim', loc(2) + span};
					m.set.axes = lim;
			end
			switch type % select plot type
				case 'contour' % contour of response to achromatic stimulus
					m.prep.prep = {'image', 'contour'};
					m.plot.funFun = @(d, m) contour(d.x, d.y, shiftdim(d.resp, 2), ...
						.2 * [-1, 1]);
					m.set.colormap = 'turbo';
				case 'mono' % monochromatic image
					m.prep.prep = {'image', 'mono'};
					m.plot.funFun = @(d, m) imagesc(d.x, d.y, shiftdim(d.resp, 2));
					lim = 1; m.set.caxis = lim * [-1, 1];
					m.set.colorbar = {'ticks', lim * [-1, 0, 1]};
					r = linspace(0, 1, 256)'; c = [r, flip(r), 0 * r]; m.set.colormap = c;
				case 'true' % true colour image
					m.prep.prep = {'image', 'true'};
					m.plot.funFun = @(d, m) image(d.x, d.y, shiftdim(d.resp, 2));
			end
		case 'spec.fx.fy' % plot power spectrum of orientation preference map
			m.tasks = 'load spec desc plot stop set';
			m.x = 'freq'; m.y = 'freq';
			m.load.name = 'Orientation map.mat';
			m.plot.funFun = @(d, m) imagesc(d.freq, d.freq, shiftdim(d.spec, 1)');
			%	xs = .5 * m.p.xs / m.p.wid;
			%	lim = xs * [-1, 1]; tick = xs * [-1, 0, 1];
			%	m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			%	m.set.colorbar = {'xLim', [0, 1], 'ticks', [0, .5, 1], ...
			%	'tickLabels', {'-90', '0', '90'}};
		case 'struct.x.y' % plot cell locations and properties
			m.tasks = 'select stage prep desc plot set';
			m.x = 'x'; m.y = 'y';
			ana = ''; % analysis type: '', 'mech'
			mark = 0; % mark location(s) on the map: 0 or 1
			size = 30; % marker size (pt): Fig. 1 uses 30
			span = 1; % zoom: 1 for none, < 1 for visual field zoom factor
			m.select.phase = d.phase(end); m.select.cycle = d.cycle(end);
			m.stage.name = {'cone', 'gangOff', 'betaEx', 'betaIn'};
			m.stage.name = 'betaEx';
			m.prep.group = 'stage'; m.prep.prep = 'scatter';
			m.plot.group = 'stage';
			m.plot.funFun = @(d, m) ...
				scatter(d.x.val, d.y.val, size, d.col.val, 'filled');
			lim = .2; % lim = m.p.wid;
			m.set.axes = setFig(m, xSym = lim, ySym = lim);
			m.set.colorbar = {};
			switch ana
				case 'mech' % add centre and surround mechanisms
					m.tasks = replace(m.tasks, 'plot', 'plot fig');
					m.fig.mech = m.p.c.e; % mechanism location
			end
			switch mark % mark a location on the map
				case 1
					m.tasks = replace(m.tasks, 'set', 'set add');
					loc = [m.p.c.a; m.p.c.b; m.p.c.c; m.p.c.d];
					m.add.funFun = @(d, m) plot(loc(:, 1), loc(:, 2), 'ok', ...
						lineWidth = 2);
			end
			if span < 1 % zoom into cell of interest
				m.set.axes = setFig(m, zoom = {loc, span});
			end
		case 'stim.x.y' % show stimulus
			m.tasks = 'stim desc plot set';
			move = 0; % show a movie: 0 or 1
			m.p.stimS = 'grating'; % m.p.stimS = 'image';
			m.plot.funFun = @(d, m) image(d.x, d.y, ...
				shiftdim(d.image(1, :, :, :, 1), 1));
			m.set.axes = setFig(m, xSym = m.p.wid, ySym = m.p.wid);
			m.set.axes = setFig(m, xSym = .2, ySym = .2); % talk
			switch move
				case 1
					m.tasks = replace(m.tasks, 'set', 'set move');
			end
			switch m.p.stimS
				case 'grating' % set stimulus parameters
					m.p.cont = 1 * [1, -1, 0]; % contrast
					m.p.dir = 0; m.p.freqS = 10; m.p.phaseS = 0;
					%	m = setFig(m, fig2 = 1); % dir. and spat. frequency: stim. 1, 2, 3
				case 'image'
					m.tasks = replace(m.tasks, 'stim', 'load stim');
					m.load.name = 'Image 8.mat';
					m.stim.image = 2;
			end
		case 'weight.x.y' % plot synaptic weights
			m.tasks = 'select stage addW prep desc plot set add pub';
			m.x = 'x'; m.y = 'y';
			loc = m.p.c.b; % location of cortical cell
			m.prep.cone = 0; % 0 for on/off colouring, 1 for cone colouring
			span = 1; % 1 for standard visual field, < 1 for visual field zoom factor
			m.select.cycle = d.cycle(end); m.select.phase = d.phase(end);
			m.stage.name = {'genOff', 'genOn'};
			m.addW.group = 'stage'; m.addW.loc = loc;
			m.prep.group = 'stage'; m.prep.prep = 'scatter'; m.prep.diam = 2000;
			m.plot.group = 'stage';
			m.plot.funFun = @(d, m) scatter(d.x.val, d.y.val, ...
				d.size.val, d.col.val, 'filled', 'markerEdgeColor', 'k');
			lim = .2; m.set.axes = setFig(m, xSym = lim, ySym = lim);
			m.add.funFun = @(d, m)scatter(loc(1), loc(2), 50, 'k', 'filled');
			if span < 1 % zoom into cell of interest
				m.set.axes = setFig(m, zoom = {loc, span});
			end
	end

	% Initiate analysis
	m = prim.setHandle(m); % set task function handles
	m = m.calStruct(m); % calculate radii, arrays, and convergence functions
	m = m.setSyn(m); % create list of all plastic synapse
	m = m.samp(m); % create spatial and temporal samples
	stream(d, m); % execute the task list

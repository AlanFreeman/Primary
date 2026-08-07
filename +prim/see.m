function out = see(m, opt) % show interim or additional analysis results

	arguments

		m struct % metadata
		opt.add double % add a point to an existing map
		opt.format string = "image" % image format: image, imagesc, or imshow
		opt.image % show an image
		opt.index (1, 2) cell % find location index
		opt.locate % find the location of a datatip on an existing plot
		opt.movie % show a movie of images
		opt.scatter % show locations and weights of cortical input
		opt.weight % input weights

	end

	% Perform the required examination
	import prim.showStim % find function
	name = string(fieldnames(opt))'; % option field names: 1 x os, string
	for nameC = name % look for key option
		switch nameC % the first option indicates the type of examination
			case "add" % add a point to an existing map *** check ***
				a = gca; % current axes
				a.NextPlot = "add"; % add to existing axes, don't replace
				plot(a, loc(1), loc(2), "ok"); % add point
			case "image" % show an image
				im = opt.image(:, :, :, 1); % image: xs x ys x cs
				f = opt.format; % image format: image, imagesc, or imshow
				figure("windowStyle", "docked"); % assume figures window
				switch f % select format
					case "image" % image in image format
						im = permute(im, [2, 1, 3]); % swap x and y
						im = flipud(im); % flip up for down
						image(im); % show
					case "imageCont" % image in contrast units
						im = permute(im, [2, 1, 3]); % swap x and y
						im = flipud(im); % flip up for down
						im = .5 * (im + 1); % convert contrast to a range of ~ [0, 1]
						imagesc(im); % show
					case "imshow" % image in imshow format
						imshow(im); % show
				end
			case "index" % find location index
				val = opt.(nameC); % values: 1 x 2, cell
				s = val{1}; % stage
				loc = val{2}; % location
				out = knnsearch(m.p.(s).loc, loc);
			case "locate" % find the cell location of a datatip on an existing plot
				a = gca; % current axes
				o = findobj(a, type = opt.locate); % object containing datatip
				cAll = [o.XData', o.YData']; % coordinates of all points in plot: ls x 2
				d = findobj(a, type = 'datatip'); % datatip object
				c = [d.X, d.Y]; % coordinates of datatip: 1 x 2
				i = knnsearch(cAll, c); % index of datatipped point
				loc = o.UserData; % locations of all cells (deg): ls x 2
				out = loc(i, :); % required cell location (deg): 1 x 2
			case "movie" % show a movie of images
				figure("windowStyle", "docked"); % assume figures window
				switch "loop" % select method
					case "loop" % loop over frames *** check ***
						im = opt.movie; % image: xs x ys x cs x ts
						ts = size(im, 4); % number of frames
						stim.dur = m.p.time / ts; % frame duration (s)
						stim.stim = im; % stimulus: xs x ys x cs x ts
						stim.x = .5 * m.p.wid * [-1, 1]; % x limits (deg): 1 x 2
						showStim(stim); % show movie
					case "movie" % use Matlab movie function

						%	Initialise
						mov = opt.movie; % image: xs x ys x cs x fs
						if class(mov) == "double" % mov has double class
							mov = min(max(mov, 0), 1); % restrict to [0, 1] for im2frame
						end
						fs = size(mov, 4); % number of frames
						f(fs) = struct("cdata", [], "colormap", []); % allocate storage
		
						%	Create movie
						for i = 1: fs % store frames
							im = mov(:, :, :, i); % current frame
							im = permute(im, [2, 1, 3]); % swap x and y
							im = flipud(im); % flip up for down
							f(i) = im2frame(im); % create movie
						end
		
						%	Show movie
						a = gca; % axes created by figure
						axis tight manual; % allow setting of axis limits
						a.XLim = [1, 167]; a.YLim = [1, 167]; % set axis limits
						movie(a, f); % show movie

				end
			case "scatter" % show the locations and weights of cortical input
				locSub = opt.scatter{1}; % locations of subcortical input: ss x 2
				locCort = opt.scatter{2}; % locations of cortical cells: cs x 2
				locInt = opt.scatter{3}; % location of cortical cell of interest: 1 x 2
				w = opt.weight; % input weights: cs x ss
				i = knnsearch(locCort, locInt); % index of cell of interest
				w = w(i, :); % weights of interest: 1 x ss
				w = 60 * w / max(w); % prepare for plot
				figure("windowStyle", "docked"); % new figure
				scatter(locSub(:, 1), locSub(:, 2), w, 'filled') % scatter plot
		end
	end
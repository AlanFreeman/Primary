function setModel

% Calculate visual system model parameters from published work

	%	Initialise
	m.folder = userpath + "/Data/Model"; % data folder
	m = model.setLit(m); % set metadata calculated from the literature
	m = model.setPers(m); % set personal preferences

	% Set metadata for stream
	switch "count.bpi" % select analysis
		case "count.bpi" % plot "bandpass index": Lee (12)
			m.tasks = "bpi plot set pub";
			m.x = "bpi"; m.y = "gen";
			m.plot.group = "stim";
			m.plot.funFun = @(d, m) stairs([d.bpi; 1], [d.gen; 0] / max(d.gen));
			m.set.axes = prim.setFig(m, "xZero", 1, "yZero", 1);
		case "ecFunc.ecAna" % plot and print func. v. ana. eccentricity
			m.tasks = "ecFunc plot set";
			m.x = "ecAna"; m.y = "ecFunc";
			ana = ""; % analysis type: "", "check", "print"
			m.ecFunc.cut = 5; % low-eccentricity cutoff for Schein data (deg)
			m.ecFunc.inter = [.8, 0]; % low-ec. intercept for McGregor data (deg)
			m.set.axes = prim.setFig(m, "xZero", 10, "yZero", 10);
			switch ana % which analysis type?
				case "check" % check readTab output
					ecA = linspace(0, 10, 101); ecF = model.readTab("ecFunc", ecA);
					m.plot.funFun = @(d, m) plot(ecA, ecF);
				case "print" % print conversion function in preparation for readTab
					m.tasks = replace(m.tasks, "plot set", "print");
					m.print.var = "ecFunc";
			end
		case "ecFunc.ecAna low" % plot func. v. ana. eccentricity: McGregor (18)
			m.tasks = "ecLow plot set";
			m.x = "ecAna"; m.y = "ecFunc";
			ana = "fit"; % analysis type: "", "fit", "show"
			m.plot.line = "source"; m.plot.arg = {"o"};
			m.fitlm.arg = {"poly3", "predictorVars", m.x, "responseVar", m.y};
			m.set.axes = prim.setFig(m, "xZero", 3, "yZero", 1);
			switch ana % which analysis type?
				case "fit" % fit a regression
					m.tasks = replace(m.tasks, "plot", "plot fitlm pred add");
				case "show" % show regression parameters
					m.tasks = replace(m.tasks, "plot set", "fitlm show");
			end
		case "ecFunc.ecAna high" % plot func. v. ana. eccentricity: Schein (18)
			m.tasks = "ecHigh plot set";
			m.x = "ecAna"; m.y = "ecFunc";
			m.set.axes = prim.setFig(m, "xZero", 15, "yZero", 15);
		case "dens.ec" % plot ganglion cell density v. func. ec.: Wässle (89)
			m.tasks = "densGangAna densGang plot set";
			m.x = "ec"; m.y = "densOn";
			ana = "print"; % analysis type: "", "check", "print"
			m.plot.funFun = @(d, m) plot(d.ec, d.densOff, d.ec, d.densOn);
			m.set.axes = {"xScale", "log", "yScale", "log"};
			switch ana % which analysis?
				case "check" % check readTab output
					e = [0: .25: .75, 1: .5: 3.5, 4: 1: 9, 10: 5: 50];
					densOff = model.readTab("densGangOff", e);
					densOn = model.readTab("densGangOn", e);
					m.plot.funFun = @(d, m) plot(e, densOff, e, densOn);
				case "print"
					m.tasks = replace(m.tasks, "plot set", "print");
					m.print.var = "densOn"; m.print.vals = 6;
			end
		case "dens.ecAna" % plot ganglion cell density v. ana. ec.: Wässle (89)
			m.tasks = "densGangAna plot set";
			m.x = "ecAna"; m.y = "dens";
			m.plot.line = "source"; m.plot.arg = "o";
			m.set.axes = {"xScale", "log", "xLim", [.3, 60], ...
				"xTick", [.3, 1, 10, 60], "yScale", "log", "yLim", [3, 3000]};
		case "rad.ec gang" % plot and print gang. cell dendritic radius v. func. ec.
			m.tasks = "radGang plot set";
			m.x = "ec"; m.y = "radius";
			ana = "check"; % analysis type
			switch ana % which analysis?
				case "check" % check readTab output
					eStart = m.p.ecTrans(1); eEnd = m.p.ecTrans(2); % start, end of trans.
					e = [0, eStart: eEnd - 1, eEnd: 5: 60]; % eccentricity (deg): 1 x es
					radius = model.readTab("radGang", e, m);
					m.plot.funFun = @(d, m) plot(e, radius);
				case "print"
					m.tasks = replace(m.tasks, "plot set", "print");
					m.print.var = "radius";
			end
		case "rad.ec hor" % plot and print horizontal cell rec. field radius v. ec.
			m.tasks = "radHor select plot set";
			m.x = "ec"; m.y = "radius";
			ana = "plot"; % analysis: "check", "fit", "plot", "print"
			m.select.type = "H2";
			m.plot.group = "type"; m.plot.arg = "o";
			%	m.set.axes = prim.setFig(m, "xZero", 70, "yZero", 1.3);
			m.set.axes = {"xLim", [0, 60], "yScale", "log", "yLim", [.03, 1.3], ...
				"clipping", "off"};
			switch ana % choose analysis
				case "check" % check readTab output
					e = linspace(0, 50, 101); r = model.readTab("radHor2", e);
					m.plot.funFun = @(d, m) plot(e, r);
				case "fit"
					m.tasks = replace(m.tasks, "set", "radHorFit add set");
				case "print"
					m.tasks = replace(m.tasks, "plot set", "radHorFit print");
					m.print.var = "radius";
			end
	end

	% Run tasks
	d = table; % empty data table
	list = stream.listTasks("model", [], m); % list task handles and grouping var.
	stream(d, m, list); % run tasks
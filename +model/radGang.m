function [d, m] = radGang(~, m) % ganglion cell dendritic radius: Watanabe (89)

	% Load diameter data and convert both eccentricity and diameter
	folder = m.folder + filesep + "Watanabe (89)"; % data folder
	d = prim.readFolder(folder); % create table from folder
	d.label = []; % remove this variable: it is not applicable
	d = sortrows(d, "eccMm"); % sort on eccentricity
	ecAna = m.p.magRet * d.eccMm; % anatomical eccentricity (deg): es x 1
	ecEmp = model.readTab("ecFunc", ecAna); % convert ana. to func. ec. (deg)
	diam = 10 .^ d.diamLog; % dendritic field diameter (um): es x 1
	r = (.5 * diam / 1000); % radius (mm): es x 1
	r = m.p.magRet * r; % radius (deg): es x 1

	%	Resample for efficient storage and set radius for private line
	eStart = m.p.ecTrans(1); eEnd = m.p.ecTrans(2); % start, end of transition
	ec = [0, eStart: eEnd - 1, eEnd: 5: 60]'; % eccentricity (deg): es x 1
	rPriv = 0 * ec; % private line radii (deg): es x 1
	
	%	Resample empirical data and merge private-line and empirical data
	modelEmp = fitlm(ecEmp, r, "poly3"); % fit empirical data
	rEmp = predict(modelEmp, ec); % fitted data
	radius = model.trans(ec, rPriv, rEmp, m); % zero radius at low eccentricity

	% Store
	d = table(ec, radius); % store in table
	d.Properties.VariableDescriptions = ...
		["Eccentricity (deg)", "Dendritic field radius (deg)"];
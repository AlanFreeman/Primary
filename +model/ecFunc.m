function [d, m] = ecFunc(~, m) % make lookup table to convert anat. to func. ec.

	%	Fit polynomial to low eccentricity data: McGregor (18)
	[dLow, m] = model.ecLow([], m); % low eccentricity data
	modLow = fitlm(dLow, "poly3", ...
		"predictorVars", "ecAna", "responseVar", "ecFunc"); % fit
	ecA = dLow.ecAna; % anatomical eccentricity (deg): es x 1
	ecF = predict(modLow, dLow); % functional eccentricity (deg)
	ecLow = [ecA, ecF]; % anat. and func. eccentricity (deg): es x 2

	%	Calculate and cut high eccentricity data: Schein, 88
	[dHigh, m] = model.ecHigh([], m); % high eccentricity data
	ecA = dHigh.ecAna; % anatomical eccentricity (deg): es x 1
	ecHigh = [ecA, dHigh.ecFunc]; % combine anat., func.: es x 2
	i = ecA >= m.ecFunc.cut; % indices of data above cut point
	ecHigh = ecHigh(i, :); % keep only those data: es x 2

	%	Fit spline to combined data
	ec = [ecLow; ecHigh]; % combine low and high eccentricity data (deg): es x 2
	inter = m.ecFunc.inter; % low eccentricity intercept
	ec = [inter; ec]; % add low eccentricity intercept (deg): es x 2
	junc = max(ecLow(:, 1)); % top end of low eccentricity data;
	xLow = linspace(inter(1), junc, 10)'; % low interpolation points
	xHigh = linspace(junc, max(ecHigh(:, 1)), 11)'; % high interp. points
	xHigh = xHigh(2: end, :); % remove overlapping point
	x = [xLow; xHigh]; % interpolation points (deg)
	y = spline(ec(:, 1), ec(:, 2), x); % fit spline
	ec = [x, y]; % interpolated data: es x 2

	%	Prepare data for storage in readTab
	ec = [[0, 0]; ec; [15, 15]; [90, 90]]; % add endpoints (deg): es x 2
	ecAna = ec(:, 1); ecFunc = ec(:, 2); % prepare to store
	d = table(ecAna, ecFunc); % store
	d.Properties.VariableDescriptions("ecAna") = "Anatomical eccentricity (deg)";
	d.Properties.VariableDescriptions("ecFunc") = "Functional eccentricity (deg)";
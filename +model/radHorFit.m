function [d, m] = radHorFit(d, m) % fit hor. cell convergence function radius

	%	Initialise
	type = d.type(1); % horizontal cell type
	if type == "H1" % H1 horizontal cell
		eTrans = [16, 20]; % transition range from low to high eccentricity (deg)
		ec = [0: 5, 6: 2: 14, 16: 19, 20: 5: 50]';
			% eccentricity for efficient storage in readTab (deg): es x 1
		mult = 1; % multiplier of high-eccentricity data
	else % H2 horizontal cell
		eTrans = [20, 30];
		ec = [0: 4: 16, 20: 28, 30: 5: 50]'; % eccentricity (deg): es x 1
		mult = .37; % H2s have smaller receptive fields, Zhang (11) Light
	end

	%	Fit the low-eccentricity data
	e = d.ec; r = d.radius; % data to fit
	modelLow = fitlm(e, r, "poly3"); % fit empirical data
	r = predict(modelLow, e); % fitted data

	%	Select and store the low-eccentricity data
	i = e < eTrans(1); % indices of low-eccentricity data
	eLow = e(i); rLow = r(i); % select
	
	%	Generate the high-eccentricity data
	eHigh = linspace(eTrans(2), ec(end), 20)'; % finely spaced eccentricities
	rHigh = exp(polyval(m.p.radHorCoef, eHigh));
		% use polynomial fitted to H1 horizontal cell receptive field radii (deg)
	rHigh = mult * rHigh; % reduce H2 radii

	%	Join the low- and high-eccentricity data
	e = [eLow; eHigh]; r = [rLow; rHigh]; % join
	r = spline(e, r, ec); % reduce the number of points and smooth
	es = length(ec); % number of eccentricities
	type = repmat(type, [es, 1]); % for storage
	d = table(type, ec, r, VariableNames = ["type", "ec", "radius"]); % store
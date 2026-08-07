function p = readTab(name, ecc, m) % interp. model parameter at specified ecc.

	switch name
		case 'densCone' % cone density
			e = [
 				 0.0000,  0.3017,  0.5465,  0.8296,  1.0674,  1.5659,  2.0644, ...
 				 2.5896,  3.0588,  3.5418,  4.0933,  4.5761,  5.0863, 10.1849, ...
				15.2140, 20.6450, 25.6239, 30.5727, 35.6190, 40.8653, 45.6091, ...
				50.7488, 55.7820, 61.0172, 65.8515, 71.0757, 76.2114
				]; % eccentricity (deg): 1 x es
			r = [
				8930.0976, 6587.4855, 4927.0627, 3730.0137, 3150.7964, 2713.1656, ...
				2262.6633, 1953.7474, 1786.4179, 1606.2170, 1477.5021, 1323.0441, ...
				1207.2007,  662.7883,  509.1457,  355.5032,  275.0238,  238.4422, ...
 				 206.2504,  175.5219,  144.7934,  122.8445,  109.6751,   95.0425, ...
  				78.9466,   80.4099,   64.3140
				]; % density (cells/deg^2): 1 x es
		case 'densGangOff' % off-centre ganglion cell density
			e = [0: .25: .75, 1: .5: 2.5, 3: 1: 9, 10: 5: 50]; % ecc. (deg): 1 x es
			r = [
				8930.0976, 6988.9209, 5242.4616, 4066.5912, 3314.9648, 2284.4658, ...
				1658.0154, 1314.5218, 1132.8928,  915.2148,  770.2689,  658.9283, ...
				 565.7395,  486.5021,  415.5119,  352.7161,  184.7380,   93.5624, ...
  				51.8562,   31.9776,   19.8953,   13.4204,    8.7238,    6.3519
				]; % density (cells/deg^2): 1 x es
		case 'densGangOn' % on-centre ganglion cell density
			e = [0: .25: .75, 1: .5: 2.5, 3: 1: 9, 10: 5: 50]; % ecc. (deg): 1 x es
			r = [
				8930.0976, 6988.9209, 5242.4616, 4066.5912, 3314.9648, 2269.7645, ...
				1615.3249, 1243.4008, 1045.7472,  844.8136,  711.0175,  608.2415, ...
				 522.2211,  449.0788,  383.5494,  325.5841,  170.5274,   86.3653, ...
  				47.8673,   29.5178,   18.3649,   12.3880,    8.0528,    5.8633
				]; % density (cells/deg^2): 1 x es
		case 'devGang' % standard deviation of ganglion cells about triangular array
			e = [0, 1, 3, 50]; % eccentricity (deg): 1 x es
			dev = m.p.kGangDev; % standard deviation (deg)
			r = [0, 0, dev, dev]; % deviation (deg): 1 x es
		case 'freqS' % spatial frequency cutoff
			% Use setPrim > freq.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [24.75, 18, 9.9, 3.6, .9]; % spatial frequency cutoff (cycles/deg)
		case 'radGang' % ganglion cell dendritic field radius
			%{
			e = [0, 1, 2, 3: 3: 60]; % eccentricity (deg): 1 x es
			r = [
				0.0010, 0.0044, 0.0079, 0.0113, 0.0208, 0.0320, 0.0448, 0.0590, ...
				0.0744, 0.0908, 0.1082, 0.1262, 0.1448, 0.1637, 0.1828, 0.2019, ...
				0.2209, 0.2395, 0.2577, 0.2751, 0.2917, 0.3073, 0.3217
				]; % radius (deg): 1 x es
			%}
			e = [0, 1, 9: 3: 60]; % eccentricity (deg): 1 x es
			r = [
				0.0000, 0.0000, 0.0346, 0.0461, 0.0594, 0.0743, 0.0904, 0.1076, ... 
				0.1257, 0.1445, 0.1636, 0.1830, 0.2024, 0.2215, 0.2402, 0.2582, ...
				0.2753, 0.2913, 0.3059, 0.3190
				]; % radius (deg): 1 x es
		case 'ratCort' % inverse cortical magnification factor
			e = linspace(0, 3000, 20) / 60; % eccentricity (deg): 1 x es
			r = [
				0.0854, 0.2914, 0.5125, 0.7544, 1.0221, 1.3174, 1.6416, 1.9958, ...
				2.3807, 2.7973, 3.2464, 3.7289, 4.2456, 4.7975, 5.3853, 6.0100, ...
				6.6725, 7.3737, 8.1145, 8.8958
				]; % % inverse cortical magnification factor (deg/mm)
		case 'sepEx' % excitatory cell separation
			% Use setPrim > sep.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			%	r = [.005, .008, .015, .026, .09]; % sep., fixed ratio (deg): 1 x es
			r = [.005, .008, .01, .016, .028]; % sep., fixed number (deg): 1 x es
		case 'sepIn' % inhibitory cell separation
			% Use setPrim > sep.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			%	r = [.012, .02, .036, .06, .21]; % sep., fixed ratio (deg): 1 x es
			r = [.012, .02, .025, .039, .07]; % sep., fixed number (deg): 1 x es
		case 'width' % visual field width required for fixed number of cones
			% Use setPrim > width.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [.15, .25, .32, .5, .9]; % visual field width (deg)
	end
	p = interp1(e, r, ecc); % interpolate radius at specified ecc. (deg)

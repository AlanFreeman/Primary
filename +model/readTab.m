function p = readTab(name, ec, m) % interp. model parameter at specified ec.

	switch name
		case "densCone" % cone density
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
		case "densGangOff" % off-centre ganglion cell density
			e = [0: .25: .75, 1: .5: 3.5, 4: 1: 9, 10: 5: 50]; % ec. (deg): 1 x es
			r = [
				8930.0976, 6988.9209, 5242.4616, 4066.5912, 3314.9648, 2771.0189, ...
				2320.8626, 2006.4490, 1807.3876, 1621.8120, 1499.2774, 1187.3535, ...
				 997.5529,  809.9393,  647.8856,  509.1454,  395.4397,  229.9994, ...
				 111.6385,  54.0890,   31.4514,   20.8822,   13.5861,     6.7343, ...
				   4.2464
				]; % density (cells/deg^2): 1 x es
		case "densGangOn" % on-centre ganglion cell density
			e = [0: .25: .75, 1: .5: 3.5, 4: 1: 9, 10: 5: 50]; % ec. (deg): 1 x es
			r = [
				8930.0976, 6988.9209, 5242.4616, 4066.5912, 3314.9648, 2771.0189, ...
				2320.8626, 2006.4490, 1807.3876, 1621.8120, 1499.2774, 1132.3455, ...
				 911.6216,  707.1032,  531.4887,  383.2807,  263.6265,  153.3329, ...
				  74.4257,   36.0593,   20.9676,   13.9215,    9.0574,    4.4895, ...
				   2.8309
				]; % density (cells/deg^2): 1 x es
		case "devGang" % standard deviation of ganglion cells about triangular array
			e = [0, m.p.ecTrans, 50]; % eccentricity (deg): 1 x es
			dev = m.p.kGangDev; % standard deviation (deg)
			r = [0, 0, dev, dev]; % deviation (deg): 1 x es
		case "ecFunc" % conversion of anatomical to functional eccentricity
			e = [
				0.0000, 0.8000,  1.0097,  1.2194,  1.4291,  1.6389,  1.8486, 2.0583, ...
				2.2680, 2.4777,  2.6874,  3.6243,  4.5611,  5.4979,  6.4347, 7.3715, ...
				8.3083, 9.2451, 10.1820, 11.1188, 12.0556, 15.0000, 90.0000
				]; % anatomical eccentricity (deg): 1 x es
			r = [
				0.0000, 0.0000, 0.0201,  0.0339,  0.0694,  0.1252,  0.1995, 0.2907, ...
				0.3970, 0.5171, 0.6485,  1.6920,  3.0758,  4.2296,  5.3142, 6.4194, ...
				7.5474, 8.6544, 9.7376, 10.7735, 11.7594, 15.0000, 90.0000
				]; % functional eccentricity (deg): 1 x es
		case "freqS" % spatial frequency cutoff
			% Use setPrim > freq.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [24.75, 18, 9.9, 3.6, .9]; % spatial frequency cutoff (cycles/deg)
		case "radGang" % ganglion cell dendritic field radius
			eStart = m.p.ecTrans(1); eEnd = m.p.ecTrans(2); % start, end of transition
			e = [0, eStart: eEnd - 1, eEnd: 5: 60]; % eccentricity (deg): 1 x es
			r = [
				0.0000, 0.0000, 0.0037, 0.0083, 0.0140, 0.0208, 0.0288, 0.0382, ...
				0.0595, 0.0849, 0.1136, 0.1445, 0.1765, 0.2088, 0.2402, 0.2697, ...
				0.2963, 0.3191
				]; % radius (deg): 1 x es
		case "radHor" % H1 horizontal cell receptive field radius
			e = [0: 5, 6: 2: 14, 16: 19, 20: 5: 50]; % eccentricity (deg): 1 x es
			r = [
				0.0442, 0.0457, 0.0474, 0.0492, 0.0512, 0.0533, 0.0556, 0.0605, ...
				0.0660, 0.0720, 0.0785, 0.0855, 0.0974, 0.1182, 0.1427, 0.1656, ...
				0.2709, 0.4357, 0.6650, 0.9313, 1.1562, 1.2296
				]; % radius (deg): 1 x es
		case "radHor2" % H2 horizontal cell receptive field radius
			e = [0: 4: 16, 20: 28, 30: 5: 50]; % eccentricity (deg): 1 x es
			r = [
				0.0760, 0.0834, 0.0907, 0.0979, 0.1050, 0.1120, 0.1134, 0.1147, ...
				0.1162, 0.1182, 0.1211, 0.1252, 0.1309, 0.1386, 0.1612, 0.2461, ...
				0.3446, 0.4278, 0.4550
				]; % radius (deg): 1 x es
		case "ratCort" % inverse cortical magnification factor
			e = linspace(0, 3000, 20) / 60; % eccentricity (deg): 1 x es
			r = [
				0.0854, 0.2914, 0.5125, 0.7544, 1.0221, 1.3174, 1.6416, 1.9958, ...
				2.3807, 2.7973, 3.2464, 3.7289, 4.2456, 4.7975, 5.3853, 6.0100, ...
				6.6725, 7.3737, 8.1145, 8.8958
				]; % % inverse cortical magnification factor (deg/mm)
		case "sepEx" % excitatory cell separation
			% Use setPrim > sep.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [.005, .008, .01, .016, .028]; % sep., fixed number (deg): 1 x es
		case "sepIn" % inhibitory cell separation
			% Use setPrim > sep.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [.012, .02, .025, .039, .07]; % sep., fixed number (deg): 1 x es
		case "width" % visual field width required for fixed number of cones
			% Use setPrim > width.ecc to set these (pragmatic) values
			e = [0, 1, 3, 10, 30]; % eccentricity (deg): 1 x es
			r = [.15, .25, .32, .5, .9]; % visual field width (deg)
	end
	p = interp1(e, r, ec); % interpolate radius at specified ec. (deg)
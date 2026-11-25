function m = gen(m)

%	Create geniculate array

	m.p.gen.loc = [m.p.gangOff.loc; m.p.gangOn.loc];
		% LGN locations (x, y) (deg): ls x 2
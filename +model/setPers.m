function m = setPers(m)

%	Set personal preferences

	%	Analysis variables
	m.x = "x"; % default independent variable
	m.y = "resp"; % default independent or dependent variable
	m.z = "resp"; % default dependent variable
	
	%	Publication variables
	m.pub.units = "centimeters"; % define unit before using it
	m.pub.axisLength = [5, 5]; % length of x-, y-axes (cm)
	m.pub.dest = "pdf"; % export destination: "copy", "eps", or "pdf"
	m.pub.lineWidth = 1.5; % width of axes and lines (cm)
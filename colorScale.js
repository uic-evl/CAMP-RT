var ColorScale = (function(){
	
	var numSegments = 20;
	
	var drawScale = function(target, colorMap, domain){
		var div = document.getElementById(target);
		var height = div.clientHeight;
		var width = div.clientWidth;
		var binHeight = height/numSegments;
		var doseInterval = (domain[1] - domain[0])/numSegments;
		var colors = [];
		var currentDose = domain[0];
		for(var i = 0; i < numSegments; i++){
			var newColor = colorMap(currentDose);
			currentDose = currentDose + doseInterval;
			colors.push(newColor);
		}
		console.log(domain);
		var svg = d3.select('#'+target).append('svg')
			.attr('height', height)
			.attr('width', width);
		svg.selectAll('rect')
			.data(colors).enter()
			.append('rect')
			.attr('x', 0)
			.attr('y', function(d,i){ return i*binHeight;})
			.attr('fill', function(d){ return d; })
			.attr('width', width)
			.attr('height', binHeight)
			.attr('opacity', .5);
	}
	
	var draw = function(){
		drawScale('doseColorScale', Controller.getDoseColor, data.getMeanDoseExtents());
		drawScale('doseErrorColorScale', Controller.getDoseErrorColor, data.getDoseErrorExtents());
	}
	
	return {
		draw: draw
	}
})();
var Controller = (function(){
	var scatterData = {};
	var bubbleData = {};
	var enableBrush = false;
	var doseColor = d3.scaleLinear()
		.domain([0, 105])
		.range(['#ffffe0','#8b0000']);
	var doseErrorColor = d3.scaleLinear()
		.domain([0, 20])
		.range(['#999999','#3d32ff']);

	return {
		
		getDoseColor: function(d){ return doseColor(d); },
		
		getDoseErrorColor: function(d){ return doseErrorColor(d); },
		
		toggleBrush: function(isEnabled){
			enableBrush = isEnabled;
		},
		
		brush: function(id){
			if(enableBrush){
				//highlight patient matches in other views on mousover
				var patientView = document.getElementById(id);
				if(patientView != null){
					var description = patientView.querySelector('.description');
					description.style.color = 'white';;
				}
				
				var bubble = d3.selectAll('#organBubble'+id);
				if(!bubble.empty()){
					bubbleData.fill = bubble.attr('fill');
					bubbleData.stroke = bubble.attr('stroke');
					bubbleData.opacity = bubble.attr('opacity');
					bubble.attr('fill', 'white')
						.attr('stroke', 'white')
						.attr('opacity', 1)
						.moveToFront();
				}
				var scatterSelection =  d3.selectAll('path').filter('#scatterDot' + id);
				if(!scatterSelection.empty()){
					scatterData.stroke = scatterSelection.attr('stroke');
					scatterSelection.attr('stroke', 'white')
						.moveToFront();
				}
			}
		},
		
		unbrush: function(id){
			if(enableBrush){
				var patientView = document.getElementById(id);
				if(patientView != null){
					var description = patientView.querySelector('.description');
					description.style.color = 'black';
				}
				if(!d3.selectAll('#organBubble'+id).empty()){
					d3.selectAll('#organBubble'+id)
						.attr('fill', bubbleData.fill)
						.attr('stroke', bubbleData.stroke)
						.attr('opacity', bubbleData.opacity);
				}
				if(!d3.selectAll('path').filter('#scatterDot' + id).empty()){
					d3.selectAll('path').filter('#scatterDot' + id)
						.attr('stroke', scatterData.stroke);
				}
				if( !d3.selectAll('.doseRect').empty()){
					d3.selectAll('.doseRect').moveToFront();
				}
			}
		},
		 
		setup: function(){
			var self = this;
			d3.selectAll('.description').on('mouseover', function(d){
				var id = this.parentNode.id;
				self.brush(id);
			}).on('mouseout', function(d){
				var id = this.parentNode.id;
				self.unbrush(id);
			});
		}
	}
})();
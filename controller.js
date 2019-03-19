var Controller = (function(){
	return {
		brush: function(id){
			//highlight patient matches in other views on mousover
			var patientView = document.getElementById(id);
			if(patientView != null){
				var description = patientView.querySelector('.description');
				description.style.color = 'white';;
			}
			d3.selectAll('circle').filter('#organBubble'+id)
				.attr('fill', 'white')
				.attr('stroke-fill', 'white')
				.attr('opacity', 1)
				.moveToFront();
			d3.selectAll('circle').filter('#scatterDot' + id)
				.attr('stroke', 'white')
				.moveToFront();
		},
		
		unbrush: function(id){
			var patientView = document.getElementById(id);
			if(patientView != null){
				var description = patientView.querySelector('.description');
				description.style.color = 'black';
			}
			d3.selectAll('circle').filter('#organBubble'+id)
				.attr('fill', 'hsl(60, 60%, 50%)')
				.attr('stroke', 'black')
				.attr('opacity', .5);
			d3.selectAll('circle').filter('#scatterDot' + id)
				.attr('stroke', 'black');
			d3.selectAll('.doseRect').moveToFront();
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
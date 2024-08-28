

function initialize_umap_view (dataset = ORIGINAL_DATASET) {

    show_view("umap");

    $UMAP.empty();

    // init global SVG and MARGIN
    UMAP_MARGIN = {top: 10, right: 60, bottom: 100, left: 0};

    // Initialize the SVG container for UMAP plot
    UMAP_SVG = d3.select("#umap-svg-container").append("svg")
        .attr("id", "umap-svg")
        .attr("width", $UMAP.width())
        .attr("height", $UMAP.height())
        .append("g")
        .attr("transform", "translate(" + UMAP_MARGIN.left + "," + UMAP_MARGIN.top + ")");


    // Initialize UMAP with the dataset
    SS = init_umap(dataset);

}


function init_umap (dataset) {

    var data_umap = dataset;
    pn = [];
    orientations = [];

    for (var j = 0; j < data_umap.length; j ++) {
    	pn[j] = data_umap[j]["Patient"];
        orientations[j] = data_umap[j]["ORIENTATION"];
    }

    var uu = [];
    var vv = [];
    for (var j = 0; j < dataset.length; j ++) {
            uu[j] = dataset[j]["u"] 
            vv[j] = dataset[j]["v"] 
        };

   	console.log('UMAP Correct!')

    // colors1 = Array(ORIGINAL_DATASET.length).fill('#00304e');

    var colors1 = orientations.map(function(orientation) {
        if (orientation == 0) {
            return '#FF0000';
        } else if (orientation == 1) {
            return '#3300FF';
        } else if (orientation == 2) {
            return '#00FF00';
        } else if (orientation == 3) {
            return '#FF00CC'
        } else {
            return '#00304e';
        }
    });
    
    // var update1 = {'marker':{color: colors1, size:10}};

    var trace = {
      x: uu,
      y: vv,
      mode: 'markers',
      type: 'scatter',
      hoverinfo: 'text',
      text: pn,
      marker: { size: 8, color: colors1}
    };

	var dataaa = [trace];

	var layout = {
      height: 700,
      xaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: false,
        showticklabels: false,
      },
      yaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: false,
        showticklabels: false,
      },
      title:'<b>UMAP Plot</b>',
      titlefont: {
        family: 'Titillium Web',
        size: 18,
        color: 'black'
      },
    };
	
    Plotly.newPlot('umap-svg-container', dataaa, layout, {showSendToCloud: true, scrollZoom: true});

    var myPlot = document.getElementById('umap-svg-container');

    // Ce qui se passe quand on clique sur un point
    myPlot.on('plotly_click', function(data){
        var pn = data.points[0].pointNumber,
        colors2 = colors1.slice();
        // colors2 = Array(ORIGINAL_DATASET.length).fill('#00304e');
        colors2[pn] = '#ffc000'; // highlight color // yellow
        
        var u1 = {'marker':{color: colors1, size:10}};
        var update2 = {'marker':{color: colors2, size:10}};
        Plotly.restyle('umap-svg-container', u1);
        Plotly.restyle('umap-svg-container', update2);
        
        enter_select_mode(data.points[0].text, true);

    });

    // Ce qui se passe quand on sélectionne un ou plusieurs points
    myPlot.on('plotly_selected', function(eventData) {
        console.log("plotly_selected triggered");
        let selectedPoints = [];
        if (eventData) {
            eventData.points.forEach(function(pt) {
                selectedPoints.push(data_umap[pt.pointIndex]["Patient"]);
            });
        }

        // On affiche les points sélectionnés dans un tableau
        const selectedPointsDiv = document.getElementById('umap-selected-points');
        if (selectedPoints.length > 0) {
            let table = '<table border = "1"><tr><th>Selected Points</th></tr>';
            selectedPoints.forEach(function(point) {
                table += `<tr><td>${point}</td></tr>`;
            });
            table += '</table>';
            selectedPointsDiv.innerHTML = table;
        } else {
            selectedPointsDiv.innerHTML = 'No points selected';
        }

    });

    return {
        // fourth: update1,
        colors: colors1
    };

}


function enter_select_umap_view (case_name) {
    exit_select_umap_view();

   var myPlott = document.getElementById('umap-svg-container');
    var datagraph = myPlott.data;

    for (var j = 0; j < datagraph[0].text.length; j ++) {
        if (datagraph[0].text[j] == case_name) {
            var colors3 = SS.colors1.slice();
            // var test_value = 1;
            // colors3 = Array(ORIGINAL_DATASET.length).fill('#00304e');
            colors3[j] = '#ffc000';
            var update3 = {'marker':{color: colors3, size:10}};
            Plotly.restyle('umap-svg-container', update3);
        } // else {
            // var test_value = 0;
        // }
    };

}



function exit_select_umap_view () {
    // Plotly.restyle('umap-svg-container', SS.fourth);
    Plotly.restyle('umap-svg-container', { 'marker.color': SS.colors1});
}


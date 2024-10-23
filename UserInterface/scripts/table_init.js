function initialize_data_table (dataset) {

	show_view("table");
	var $table = $("#result-table");
	var $table_chart = $("#result-table_chart");

	generate_table(dataset, $table);
	generate_config(dataset);                  // for Cols

	TABLE = $table.DataTable(DATA_TABLE_CONFIG);

	TABLE_CHART = $table_chart.DataTable(DATA_TABLE_CONFIG_CHART);

	init_visibility();
	// init_editability();
	init_button_style();

	CURRENT_SORT_ATTRIBUTE = ORIGINAL_FEATURE_LIST[TABLE.order()[0][0]];

	$table.find("tbody").on("click", 'td', function () {

		if ($(TABLE.column($(this).index() + ":visIdx").header()).text().trim() != "comments") {
			var case_name = $(this).parent().find("td:first-child").text();
			enter_select_mode(case_name, true);
		} else {
			$("tr.selected").removeClass("selected");
		}

	});

	$(".dataTables_scrollHeadInner > table > thead > tr > th").on("click", function () {
		data_sorting($(this).text(), (TABLE.order()[0][1] == 'desc'));
		update_views();
	});
}


function generate_table(dataset, table) {
    // Step 1: Create the table header with Image, New Tag, and Editable Tag Clone columns
    var thead_content = "<tr>";
    thead_content += "<th>Image</th><th>Tag</th><th>QC_Tag</th>"; // First columns: Image, New Tag, and Editable Clone
    ORIGINAL_FEATURE_LIST.slice(1).forEach(function(d) {  // Skip the original first column
        thead_content += "<th>" + d + "</th>";
    });
    thead_content += "</tr>";

    // Step 2: Create the table body
    var tbody_content = "";
    for (var i = 0; i < dataset.length; i++) {
        tbody_content += "<tr>";
        
        // Get the value of the first column (assuming it contains the image name)
        var imageName = dataset[i][ORIGINAL_FEATURE_LIST[0]];
        
        // Split the image name based on underscores to extract the new tag
        var imageParts = imageName.split('_');
        var newTag = imageParts.slice(4, 5).join('_') || 'N/A'; // 5th part for New Tag, or 'N/A' if missing
        
        // Add the Image column (unchanged)
        tbody_content += "<td>" + imageName + "</td>"; // Image (unchanged)
        
        // Add the New Tag column
        tbody_content += "<td>" + newTag + "</td>"; // New Tag
        
        // Add the Editable Clone of the New Tag column
        tbody_content += "<td contenteditable='true'>" + newTag + "</td>"; // Editable Clone of New Tag

        // Add the remaining columns (excluding the first one which we already processed)
        for (var j = 1; j < ORIGINAL_FEATURE_LIST.length; j++) {
            var cellContent = dataset[i][ORIGINAL_FEATURE_LIST[j]];
            if (typeof cellContent === 'number') {
                if (Math.abs(cellContent) >= 1e5) {
                    cellContent = cellContent.toExponential(2);
                } else {
                    cellContent = cellContent.toFixed(2);
                }
            } else if (cellContent === undefined || cellContent === null || cellContent === '') {
                cellContent = 'N/A';
            }

            tbody_content += "<td>" + cellContent + "</td>";
        }
        tbody_content += "</tr>";
    }

    // Step 3: Populate the table
    table.children("thead").empty().html(thead_content);
    table.children("tbody").empty().html(tbody_content);

    // Step 4: Add event listener to save changes when editing the cloned column
    table.on('blur', 'td.editable', function () {
        var newValue = $(this).text();
        var rowIndex = $(this).closest('tr').index();  // Get the row index
        var columnIndex = $(this).index();  // Get the column index (second column, index = 1)

        // Save the change to the dataset
        dataset[rowIndex][ORIGINAL_FEATURE_LIST[0]] = newValue;  // Save the edited value
        console.log("Row " + rowIndex + ", Column " + columnIndex + " updated to: " + newValue);
    });

    // Initialize the table with DataTables features (pagination, sorting, etc.)
    if ($.fn.DataTable.isDataTable(table)) {
        table.DataTable().clear().destroy(); // Destroy existing table if it's already initialized
    }
    TABLE = table.DataTable(DATA_TABLE_CONFIG);
}



function generate_config (dataset) {

	// 1. named column
	// 2. customized colvis

	var colvis_action = function (e, dt, node, config) {
		var column_name = node[0].text;
		if (this.active()) {
			// update the table column
			this.active(false);
			TABLE.column(column_name + ":name").visible(false);
			
			CURRENT_HIDDEN_COLUMNS.push(column_name);
			
			// update parallel coordinate -> delete from CURRENT_PARAC_ATTRIBUTES
			CURRENT_PARAC_ATTRIBUTES = generate_current_parac_attributes();
			update_chart_view("parallel_coordinate", CURRENT_MULTI_SELECTED);

		} else {
			// update the table column
			this.active(true);
			TABLE.column(column_name + ":name").visible(true);
			
			var index = CURRENT_HIDDEN_COLUMNS.indexOf(column_name);
			if (index > -1) {
				CURRENT_HIDDEN_COLUMNS.splice(index, 1);
			} else {
				console.log("[DEBUG] " + column_name + " is not in CURRENT_HIDDEN_COLUMNS.")
			}

			// update parallel coordinate
			CURRENT_PARAC_ATTRIBUTES = generate_current_parac_attributes();
			update_chart_view("parallel_coordinate", CURRENT_MULTI_SELECTED);

		}
	};

	DATA_TABLE_CONFIG["columns"] = [];
	var colvis_buttons_config = []; // customized colvis buttons list (every header) 

	ORIGINAL_FEATURE_LIST.forEach(function (header) {
		DATA_TABLE_CONFIG["columns"].push({
			name: header
		});
		colvis_buttons_config.push({
			text: header,
			// display: none,
			className: DEFAULT_HIDDEN_COLUMNS.indexOf(header) == -1 ? 'active' : null,
			action: colvis_action
		});

	});

	var colvis_config = {
		extend: 'collection',
		text: 'Metrics',
		buttons: colvis_buttons_config,
		fade: 500
	};

	DATA_TABLE_CONFIG["buttons"].push(colvis_config);
}


function init_visibility () {
	DEFAULT_HIDDEN_COLUMNS.forEach(function (hidden_header) {
		TABLE.column(hidden_header + ":name").visible(false);
	});
}

function init_button_style() {
  // Select the button container and change its class to vertical button group
  $(".table-control > div.dt-buttons").removeClass("btn-group").addClass("btn-group-vertical");
  
  // Select the buttons, change their class, and add some additional styling
  $(".table-control > div.dt-buttons > button").removeClass("btn-secondary").addClass("btn-outline-secondary").css({
    "margin-bottom": "5px",
    "color": "red"
  });
}


function select_row_in_table (case_name, from_table) {
	if (from_table) return;

	var offset = 0;

	TABLE.$("tr.selected").removeClass("selected");
	var target_index = TABLE.row(function(idx, data, node) {
		if (data[0] == case_name) {
			return true;
		} else {
			return false;
		}
	}).select().index();

	TABLE.row(target_index + offset).scrollTo();
}


function update_multi_selected_table_view (case_names) {
	TABLE.clear();
	TABLE.rows.add(CURRENT_MULTI_SELECTED.map(function(d) {return Object.values(d);})).draw();
}


function data_sorting (keyword, desc=false) {
	var compare = function (a, b) {
		if (a[keyword] < b[keyword]) {
			if (desc) {
				return 1;
			} else {
				return -1;
			}
		} else if (a[keyword] > b[keyword]) {
			if (desc) {
				return -1;
			} else {
				return 1;
			}
		} else {
			return 0;
		}
	}

	CURRENT_SORT_ATTRIBUTE = keyword;
	ORIGINAL_DATASET.sort(compare);
	// ORIGINAL_CASE_LIST = ORIGINAL_DATASET.map(function (d) {return d["Image"];});
	CURRENT_MULTI_SELECTED.sort(compare);
	CURRENT_CASE_LIST = CURRENT_MULTI_SELECTED.map(function (d) {return d["Image"];});
}


$("#save-button").on("click", function() {
    // Créer un tableau pour stocker les données mises à jour
    var updatedData = [];

    // Parcourir toutes les lignes du tableau, y compris celles qui ne sont pas visibles
    TABLE.rows().every(function() {
        var rowData = {};
        var rowNodes = this.nodes().to$(); // Récupérer l'élément jQuery pour la ligne

        rowNodes.find('td').each(function(index) {
            // Nom de la colonne correspondant à l'index
            var columnName = ORIGINAL_FEATURE_LIST[index];
            var cellValue = $(this).text(); // Récupère la valeur de la cellule
            rowData[columnName] = cellValue;
        });

        updatedData.push(rowData); // Ajoute les données de la ligne au tableau
    });

    // Mettre à jour le dataset original avec les nouvelles données
    ORIGINAL_DATASET = updatedData;

    console.log("Données sauvegardées : ", ORIGINAL_DATASET);
    alert("Modifications sauvegardées !");
});




// Fonction pour convertir le tableau visible avec les données modifiées en CSV
function tableToCSV(table) {
    var csv = [];
    
    // Récupérer les en-têtes visibles du tableau
    var headers = [];
    $(table).find('thead th').each(function() {
        headers.push($(this).text().trim());  // En-têtes visibles
    });
    csv.push(headers.join(','));  // Ajouter les en-têtes au CSV

    // Récupérer les données modifiées visibles du tableau
    $(table).find('tbody tr').each(function() {
        var row = [];
        $(this).find('td').each(function() {
            var cellValue = $(this).text().trim();  // Récupérer la valeur modifiée dans chaque cellule
            row.push(cellValue);
        });
        csv.push(row.join(','));  // Ajouter la ligne au CSV
    });

    return csv.join('\r\n');
}

// Gestionnaire d'événement pour exporter le tableau visible et modifié en CSV
$("#export-csv-button").on("click", function() {
    var csv = tableToCSV("#result-table");
    var dataStr = "data:text/csv;charset=utf-8," + encodeURIComponent(csv);
    var downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href", dataStr);
    downloadAnchorNode.setAttribute("download", "table_data.csv");
    document.body.appendChild(downloadAnchorNode); // nécessaire pour Firefox
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
});


// Appliquer des styles pour éviter l'extension lors de l'édition
$('#result-table').on('focus', 'td.editable', function() {
    $(this).css({
        'min-width': $(this).width(), // Fixe la largeur courante
        'max-width': $(this).width(),
        'overflow': 'hidden',
        'text-overflow': 'ellipsis',
        'white-space': 'nowrap'
    });
});
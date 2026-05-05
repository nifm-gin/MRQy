// Main function to initialize the table
function initialize_data_table (dataset) {

	show_view("table");     // switch to the table view in the user interface
	var $table = $("#result-table");    // select the main table element
	var $table_chart = $("#result-table_chart");    // select the chart version of the table

    // Generate the table structure and populate it with the dataset
	generate_table(dataset, $table);
    // Generate configuration for the table columns and other settings
	generate_config(dataset);

	// Initialize th DataTable with specific configurations
    TABLE = $table.DataTable(DATA_TABLE_CONFIG);
	TABLE_CHART = $table_chart.DataTable(DATA_TABLE_CONFIG_CHART);

    // Initialize column visibility, button styles, etc ...
	init_visibility();
	// Initialize editing capabilities if required
	init_button_style();

    // Set the current sort attribute based on the tables's initial order
	CURRENT_SORT_ATTRIBUTE = ORIGINAL_FEATURE_LIST[TABLE.order()[0][0]];

    // add click event listeners for table cells
	$table.find("tbody").on("click", 'td', function () {
        // If the clicked column is not "comments", enter selection mode
		if ($(TABLE.column($(this).index() + ":visIdx").header()).text().trim() != "comments") {
			var case_name = $(this).parent().find("td:first-child").text(); // get the case name
			enter_select_mode(case_name, true); // enter select mode for the row
		} else {
			$("tr.selected").removeClass("selected"); // deselect row if "comments" is clicked
		}
	});

    $table.on('click', '.visual-qc-cell', function () {
        var $cell = $(this);
        var dropdownHTML = `
            <select class="visual-qc-dropdown">
                <option value="OK">OK</option>
                <option value="WrongTag">WrongTag</option>
                <option value="WrongFOV">WrongFOV</option>
                <option value="WrongOrgan">WrongOrgan</option>
                <option value="WrongMask">WrongMask</option>
            </select>
        `;
    
        if (!$cell.find('.visual-qc-dropdown').length) {
            var currentValue = $cell.text().trim();
            $cell.html(dropdownHTML);
            $cell.find('select').val(currentValue).focus();
            //$cell.html(dropdownHTML);
            //$cell.find('select').focus().val($cell.text().trim());
            
            $cell.find('select').on('change blur', function () {
                var selectedValue = $(this).val();
                $cell.html(selectedValue); // Update cell with the selected value
            });
        }
    });

    // add sorting functionality when clicking on column headers
	$(".dataTables_scrollHeadInner > table > thead > tr > th").on("click", function () {
		data_sorting($(this).text(), (TABLE.order()[0][1] == 'desc')); // sort by columns
		update_views(); // update any dependant views
	});

    // Keyboard navigation: ArrowUp / ArrowDown to move between rows
    $(document).off("keydown.tableNav").on("keydown.tableNav", function (e) {
        if (e.key !== "ArrowUp" && e.key !== "ArrowDown") return;
        if (!CURRENT_SELECTED) return;

        e.preventDefault(); // prevent page scrolling

        // Get all currently visible & filtered rows in display order
        var rowIndexes = TABLE.rows({ order: 'current', search: 'applied' }).indexes().toArray();

        // Find the index of the currently selected row among visible rows
        var currentPos = -1;
        for (var i = 0; i < rowIndexes.length; i++) {
            var rowData = TABLE.row(rowIndexes[i]).data();
            var caseName = $(TABLE.row(rowIndexes[i]).node()).find("td:first-child").text();
            if (caseName === CURRENT_SELECTED) {
                currentPos = i;
                break;
            }
        }

        if (currentPos === -1) return; // current selection not found

        // Calculate next position
        var nextPos;
        if (e.key === "ArrowDown") {
            nextPos = Math.min(currentPos + 1, rowIndexes.length - 1);
        } else {
            nextPos = Math.max(currentPos - 1, 0);
        }

        if (nextPos === currentPos) return; // already at the boundary

        // Get the case name of the target row
        var nextCaseName = $(TABLE.row(rowIndexes[nextPos]).node()).find("td:first-child").text();

        // Trigger selection — same as a mouse click, synchronizes all views
        enter_select_mode(nextCaseName, true);

        // Scroll the DataTable body to make the selected row visible
        var $nextRow = $(TABLE.row(rowIndexes[nextPos]).node());
        var $scrollBody = $($table.closest(".dataTables_scrollBody"));
        if ($scrollBody.length) {
            var rowOffsetTop = $nextRow[0].offsetTop;
            var scrollBodyHeight = $scrollBody.height();
            var currentScrollTop = $scrollBody.scrollTop();
            // Only scroll if the row is outside the visible area
            if (rowOffsetTop < currentScrollTop || rowOffsetTop + $nextRow.outerHeight() > currentScrollTop + scrollBodyHeight) {
                $scrollBody.stop(true).animate({ scrollTop: rowOffsetTop - scrollBodyHeight / 2 }, 80);
            }
        }
    });
}
 
 
// Function to generate the table with headers and rows
function generate_table(dataset, table) {
    // Add a unique identifier to each line if not already done
    dataset.forEach((row, index) => {
        if (!row.id) {
            row.id = `row-${index}`; // Creates a unique identifier
        }
    });

    // Step 1: Create the table header with conditional columns
    var thead_content = "<tr>";
    thead_content += "<th>Image</th>";
    
    // Check if "Tag" and "QC_Tag" columns already exist in ORIGINAL_FEATURE_LIST
    var tagExists = ORIGINAL_FEATURE_LIST.includes("Tag");
    var qcTagExists = ORIGINAL_FEATURE_LIST.includes("QC_Tag");
    var visualQCExists = ORIGINAL_FEATURE_LIST.includes("Visual_QC");

    // Add "Tag" column if it doesn't already exist
    if (!tagExists) {
        thead_content += "<th>Tag</th>";
    }

    // Add "QC_Tag" column if it doesn't already exist
    if (!qcTagExists) {
        thead_content += "<th>QC_Tag</th>";
    }

    if (!visualQCExists) {
        thead_content += "<th>Visual_QC</th>";
    }
    
    // Add remaining columns from ORIGINAL_FEATURE_LIST (excluding "Image" which is already here)
    ORIGINAL_FEATURE_LIST.slice(1).forEach(function(d) {  // Skip the original first column
        thead_content += "<th>" + d + "</th>";
    });
    thead_content += "</tr>";

    // Step 2: Create the table body
    var tbody_content = "";
    for (var i = 0; i < dataset.length; i++) {
        tbody_content += `<tr data-id="${dataset[i].id}">`; // Associate the ID to each line
        
        // Get the value of the first column (assuming it contains the image name)
        var imageName = dataset[i][ORIGINAL_FEATURE_LIST[0]];
        
        // Detect tag from image name: sort KNOWN_TAGS by length descending so the most
        // specific (longest) match wins. E.g. "T2_FLAIR" is matched before "T2" or "FLAIR".
        var sortedTags = KNOWN_TAGS.slice().sort(function(a, b) { return b.length - a.length; });
        var newTag = 'N/A';
        for (var t = 0; t < sortedTags.length; t++) {
            // Match the tag as a whole word (bounded by separators or start/end of string)
            var tagRegex = new RegExp('(^|[_\\-\\.\\s])' + sortedTags[t].replace(/[-_]/g, '[_\\-]') + '([_\\-\\.\\s]|$)', 'i');
            if (tagRegex.test(imageName)) {
                newTag = sortedTags[t];
                break;
            }
        }
        
        // Add the Image column
        tbody_content += "<td>" + imageName + "</td>"; // Image
        
        // Add the Tag column
        if (!tagExists) {
            tbody_content += "<td>" + newTag + "</td>";  // Tag
        }        
        // Add the Editable QC_Tag column
        if (!qcTagExists) {
            tbody_content += "<td contenteditable='true'>" + newTag + "</td>";  // Editable QC_Tag
        }

        // Add the new Visual_QC column
        if (!visualQCExists) {
            tbody_content += "<td contenteditable='true' class='visual-qc-cell'>OK</td>"; // Default value is OK
        }

        // Add the remaining columns (excluding the first one which we already processed)
        for (var j = 1; j < ORIGINAL_FEATURE_LIST.length; j++) {
            var cellContent = dataset[i][ORIGINAL_FEATURE_LIST[j]];
            if (typeof cellContent === 'number') {
                // Format numbers with exponential notation or fixed decimal points
                if (Math.abs(cellContent) >= 1e5) {
                    cellContent = cellContent.toExponential(2);
                } else {
                    cellContent = cellContent.toFixed(2);
                }
            } else if (cellContent === undefined || cellContent === null || cellContent === '') {
                cellContent = 'N/A'; // Handle empty values
            }
            tbody_content += "<td contenteditable='true'>" + cellContent + "</td>";
        }
        tbody_content += "</tr>";
    }

    // Step 3: Populate the table
    table.children("thead").empty().html(thead_content); // add headers
    table.children("tbody").empty().html(tbody_content); // add rows

    // Step 4: Add event listener to save changes when editing the QC_tag column
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


// Function to configure column visibility and other settings for the DataTable
function generate_config (dataset) {

	var colvis_action = function (e, dt, node, config) {
        // Define the action for toggling column visibility
		var column_name = node[0].text; // Get the column name from the button text
		if (this.active()) {
			// If the column is currently visible, hide it
			this.active(false);
			TABLE.column(column_name + ":name").visible(false);
			// Add the column to the hidden columns list
			CURRENT_HIDDEN_COLUMNS.push(column_name);
			
			// update parallel coordinate -> delete from CURRENT_PARAC_ATTRIBUTES
			CURRENT_PARAC_ATTRIBUTES = generate_current_parac_attributes();
			update_chart_view("parallel_coordinate", CURRENT_MULTI_SELECTED);

		} else {
			// update the table column
			this.active(true);
			TABLE.column(column_name + ":name").visible(true);
			// Remove the column from the hidden columns list
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

    // Initialize columns and visibility settings
	DATA_TABLE_CONFIG["columns"] = [];
	var colvis_buttons_config = []; // customized colvis buttons list (every header) 

    // Configure each column based on the original feature list
	ORIGINAL_FEATURE_LIST.forEach(function (header) {
		DATA_TABLE_CONFIG["columns"].push({
			name: header // Add column names
		});
		colvis_buttons_config.push({
			text: header, // Display header as button text
			// display: none,
			className: DEFAULT_HIDDEN_COLUMNS.indexOf(header) == -1 ? 'active' : null,
			action: colvis_action // Attach the visibility toggle action
		});
	});

    // Configure the column visibility dropdown menu
	var colvis_config = {
		extend: 'collection',
		text: 'Metrics', // Dropdown menu title
		buttons: colvis_buttons_config, // Add the custom buttons
		fade: 500 // Fade effect for the dropdown
	};

    // Add the column visibility configuration to the DataTable buttons
	DATA_TABLE_CONFIG["buttons"].push(colvis_config);
}


// Function to initialize visibility of columns
function init_visibility () {
    // Hide default hidden columns
	DEFAULT_HIDDEN_COLUMNS.forEach(function (hidden_header) {
		TABLE.column(hidden_header + ":name").visible(false);
	});
}


// Function to style control buttons for a better user interface
function init_button_style() {
    // Select the button container and change its class to vertical button group
    $(".table-control > div.dt-buttons").removeClass("btn-group").addClass("btn-group-vertical");

    // Select the buttons, change their class, and add some additional styling
    $(".table-control > div.dt-buttons > button").removeClass("btn-secondary").addClass("btn-outline-secondary").css({
        "margin-bottom": "5px",
        "color": "red"
    });
}
 
 
// Function to select a row in the table by case name
function select_row_in_table (case_name, from_table) {
	if (from_table) return; // Do nothing if the selection is from another table

	var offset = 0;

	TABLE.$("tr.selected").removeClass("selected"); // Deselect currently selected rows
	var target_index = TABLE.row(function(idx, data, node) {
		if (data[0] == case_name) { // Find the row with the matching case name
			return true;
		} else {
			return false;
		}
	}).select().index(); // Select and get the index of the row

	TABLE.row(target_index + offset).scrollTo(); // Scroll to the selected row
}


// Function to update the table view with multiple selected rows
function update_multi_selected_table_view (case_names) {
	TABLE.clear(); // Clear the current table
	TABLE.rows.add(CURRENT_MULTI_SELECTED.map(function(d) {return Object.values(d);})).draw();
}


// Function to sort the dataset based on a specific column
function data_sorting (keyword, desc=false) {
    // Comparison function for sorting
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

	CURRENT_SORT_ATTRIBUTE = keyword; // Update the current sort attribute
	ORIGINAL_DATASET.sort(compare); // Sort the original dataset
	CURRENT_MULTI_SELECTED.sort(compare); // Sort the current selection
	CURRENT_CASE_LIST = CURRENT_MULTI_SELECTED.map(function (d) {return d["Image"];}); // Update the case list
}

// Trouver qqch pour que la colonne id ne soit pas imprimee dans le tableau exporte en tsv

// Event listener to save changes made in the table
$("#save-button").on("click", function() {
    TABLE.rows().every(function() {
        var visibleRow = $(this.node()); // Ligne visible
        var rowId = visibleRow.data("id"); // Récupérez l'ID unique
        var originalRow = ORIGINAL_DATASET.find(row => row.id === rowId); // Trouvez la ligne correspondante

        var rowIndex = this.index();
        var rowData = this.data();
        
        // Get the current values in the visible table for Tag and QC_Tag
        if (originalRow) {
            // Mettez à jour les champs dans ORIGINAL_DATASET
            originalRow.Tag = visibleRow.find("td:eq(1)").text().trim(); // Adaptez selon l'index de la colonne
            originalRow.QC_Tag = visibleRow.find("td:eq(2)").text().trim();
            originalRow.Visual_QC = visibleRow.find("td.visual-qc-cell").text().trim();
        }
    });

    console.log("ORIGINAL_DATASET updated with edited values.");
    console.log(ORIGINAL_DATASET)
    alert("Changes saved!");
});


// Function to convert the visible table into a TSV (Tab-Separated Values) format
function tableToCSV() {
    var tsv = [];

    // Forcer l'ordre des colonnes pour l'export
    var headers = ["Image", "Tag", "QC_Tag", "Visual_QC"];
    
    // Ajouter les colonnes restantes dans l'ordre d'origine
    ORIGINAL_FEATURE_LIST.forEach(header => {
        if (!headers.includes(header)) {
            headers.push(header);
        }
    });

    // Ajouter les en-têtes dans l'ordre défini
    tsv.push(headers.join('\t'));

    // Parcourir les lignes et organiser les données selon l'ordre des en-têtes
    ORIGINAL_DATASET.forEach(row => {
        var line = headers.map(header => {
            var value = row[header];
            // Gérer les valeurs manquantes ou nulles
            return value === undefined || value === null ? '' : value;
        });
        tsv.push(line.join('\t'));
    });

    return FILE_HEADER + tsv.join('\r\n');  // Return complete TSV data as a string
}


// Event listener to export the table as a TSV file
$("#export-csv-button").on("click", function() {
    var tsv = tableToCSV("#result-table"); // Convert the table to TSV format
    var dataStr = "data:text/tsv;charset=utf-8," + encodeURIComponent(tsv); // Encode the TSV
    var downloadAnchorNode = document.createElement('a'); // Create a download link
    downloadAnchorNode.setAttribute("href", dataStr); 
    downloadAnchorNode.setAttribute("download", "table_data.tsv");
    document.body.appendChild(downloadAnchorNode); // nécessaire pour Firefox
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
});


// Prevent table cells from expanding during editing
$('#result-table').on('focus', 'td.editable', function() {
    $(this).css({
        'min-width': $(this).width(), // Fix current width
        'max-width': $(this).width(),
        'overflow': 'hidden',
        'text-overflow': 'ellipsis',
        'white-space': 'nowrap'
    });
});
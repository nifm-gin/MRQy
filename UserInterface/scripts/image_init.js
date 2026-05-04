function initialize_image_view (case_list) {
	// Switch to the "image" view in the user interface.
	show_view("image");
	// Select the container for displaying images and clear its content.
	var $div = $("#overview-gallery");
	$div.empty();
	// Initialize the current case list with the original case list.
	CURRENT_CASE_LIST = ORIGINAL_CASE_LIST;
	// Add a click event listener to each image. When an image is clicked, it triggers the select mode using a portion of its source path.
	$div.children("div").children("img").click(function(){
		src_list = this.src.split('/');
		enter_select_mode(src_list[src_list.length-2]);
	});
}


// Dynamically adjust the height of the "image" view based on the window height and the visible heights of other sections in the interface.
function update_image_view_height () {
	$("#image-view").outerHeight(
			$(window).height() - 
			$("header").outerHeight(includeMargin=true) - 
			$("#table-view").outerHeight(includeMargin=true) - 
			$("#table_meas-view").outerHeight(includeMargin=true) - 
			$("#chart-view").outerHeight(includeMargin=true) 
		);
}


function enter_select_image_view (dir) {
	// Hide the image gallery and the image selection button
	$("#overview-gallery").css("display", "none");
	$("#img-select-button").css("display", "none");
	// Show the button to exit the image selection view
	$("#exit-image-select-view-btn").css("display", "block");

	// Clear the contents of the selection containers
	$("#select-candidate-container > *").remove();
	$("#select-image-container > *").remove();
	// Display the image selection view
	$("#select-image-view").css("display", "flex");

	// Prepare the container to display images related to "dir"
	var $div = $("#select-image-container");

	// Regular expression for processing image names
	var re = /\s*(?:',\s'|$)\s*|\['|'\]/;
	
	// Iterate through the dataset to find images matching the participant "dir"
	for (var j = 0; j < ORIGINAL_DATASET.length; j ++) {
		if (participant_names[j] == dir){
			console.log(dir);	// Debug: Log the selected directory
			// Sort and display the image names in the container
			var image_names_new = image_names[j].split(re);
			var collator = new Intl.Collator(undefined, {numeric: true, sensitivity: 'base'});
			var myArray = image_names_new;
			const image_names_new2 = myArray.reverse(myArray.sort(collator.compare));
			for (var i = 0; i < ORIGINAL_DATASET[j]["NUM"]; i++) {
				$div.append("<img id='exibit-img' src='" + generate_img_src(dir, CURRENT_IMAGE_TYPE,image_names_new2[i])+ "'/>");
			}
		}
	}

	$div = $("#select-candidate-container");

	// Add an event listener to show a detailed view on double-clicking an image
	$("#select-candidate-container > div > img").dblclick(function(){
		enter_detail_image_view($(this).attr("file_name"), $(this).attr("img_type"), this.src);
	});
	// Add an event listener to change the main displayed image on single click
	$("#select-candidate-container > div > img").click(function(){
		$("#exibit-img").attr("src", this.src)
						.attr("img_type", $(this).attr("img_type"));
	});
	// Enable entering the detailed image view by clicking the main image
	$("#exibit-img").click(function(){
		enter_detail_image_view($(this).attr("file_name"), $(this).attr("img_type"), this.src);
	});
}


function exit_select_image_view () {
	// Reset and hide the selection containers
	$("#select-candidate-container > *").remove();
	$("#select-image-container > *").remove();
	$("#select-image-view").css("display", "none");
	$("#exit-image-select-view-btn").css("display", "none");
	// Re-show the image gallery and the image selection button
	$("#overview-gallery").css("display", "flex");
	$("#img-select-button").css("display", "block");
}


function update_multi_selected_image_view (file_names) {
	// Iterate through all cases and hide or show images based on their presence in the list of selected file names
	ORIGINAL_CASE_LIST.forEach(function (d) {
		if (file_names.indexOf(d) == -1) {
			$("#" + ORIGINAL_CASE_DICT[d]["dom_id"]).css("display", "none");
		} else {
			$("#" + ORIGINAL_CASE_DICT[d]["dom_id"]).css("display", "flex");
		}
	});
}


function calculate_height ($div) {
	// Calculate the optimal height for thumbnails within the given container
	var num_thumbs = DEFAULT_IMAGE_EXTENSIONS.length;
	var max_width = Math.floor($div.width() / Math.ceil(num_thumbs / 2)) - 5;
	var cor_height = Math.floor(max_width / $("#exibit-img").width() * $("#exibit-img").height());
	var max_height = Math.floor($div.height() / 2) - 20;

	return Math.min(max_height, cor_height);	// Return the optimal height
}


function generate_img_src (file_name, img_type_index, image_name) {
	// Generate the full file path for an image based on its name, type, and extension
	var image_extension = DEFAULT_IMAGE_EXTENSIONS[img_type_index];
	//return [DATA_PATH + file_name + "/" +  image_name];
	//return imageRoot + file_name + "/" + image_name;
	return window.imageRoot + file_name + "/" + image_name;



}

function enter_detail_image_view (file_name, img_type, src) {
	// Update the information and styling to display an image in detail
	$("#detail-image-name > span").text(file_name);
	$("#overlay-image > figure").css("width", "auto")
		.css("background-image", "url(" + src + ")");
	$("#overlay-image > figure > img").attr("src", src);
	$("#overlay-container").css("pointer-events", "all")
		.css("opacity", 1);
		// Adjust dimensions to ensure the image fits properly within the figure
	var figure_height = $("#overlay-image > figure").height(),
		figure_width = $("#overlay-image > figure").width(),
		img_height = $("#overlay-image > figure > img").height(),
		img_width = $("#overlay-image > figure > img").width();
	if (figure_height < img_height) {
		$("#overlay-image > figure").width(img_width * (figure_height / img_height));
	}
}

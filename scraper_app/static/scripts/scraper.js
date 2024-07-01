// This event listener waits until the DOM is fully loaded before running the enclosed code
document.addEventListener("DOMContentLoaded", function(){
    // Get the select element
    var select = document.getElementById('select');

    // Add event listeners to the various elements
    document.querySelector('#select').addEventListener('click', selectFunction); // Handle changes to the select element
    document.querySelector('#search_field').addEventListener('input', searchSuggestions); // Handle input changes in the search field
    document.querySelector('#proccessJsonFile').addEventListener('click', proccesJsonFile); // Handle JSON file processing
    document.querySelector('#generateJsonFile').addEventListener('click', generateJsonFile); // Handle JSON file generation
    document.querySelector('#reset_db').addEventListener('click', reset_db); // Handle database reset
    document.querySelector('#load_new_substances').addEventListener('click', load_new_substances); // Handle loading new substances
    document.querySelector('#delete_search_results').addEventListener('click', delete_search_results); // Handle deletion of search results

    // Expand the select dropdown when mouse enters and collapse it when mouse leaves
    select.addEventListener('mouseenter', function() {
        this.size = this.options.length; // Show all options
    });
    
    select.addEventListener('mouseleave', function() {
        this.size = 1; // Show only one option at a time
    });
});

var resetInterval; // Variable to store the interval for resetting the database

// Function to handle changes in the select element
function selectFunction(){
    let select = document.querySelector('#select');
    let value = select.value;
    if (value == "molecular_mass"){
        let abweichung = document.querySelector('#abweichung');
        abweichung.style["display"] = "block"; // Show the abweichung element if the selected value is "molecular_mass"
    }
    if (value != "molecular_mass"){
        let abweichung = document.querySelector('#abweichung');
        abweichung.style["display"] = "none"; // Hide the abweichung element otherwise
    }
}

// Function to show the progress of resetting the database
function show_reset_progress(resetInterval){
    $.ajax({
        type: 'GET',
        url: "/webscraper/request_how_many_json_file", // Request progress from the server
        success: function(response) {
            $(".progressbar").width(response.progress +"%"); // Update the progress bar width
            console.log(response.progress)
            if (response.progress >= 95){
                $(".progressbar").width("0%"); 
                clearInterval(resetInterval); // Clear the interval once the progress is nearly complete
                document.querySelector('.progressbar').style['display'] = "none";
                document.querySelector('.progress_message').style["display"]= "none"; 
                document.querySelector('.show_witz').style["display"]= "none";
                try{  
                    var answer = document.querySelector('.answer');
                    answer.style["display"] = "flex";
                    new_substances_loaded("Substances are scraped, please wait until they are stored in the Database");
                } catch {
                    // Handle errors if necessary
                }
                changeBackground(); // Change the background
                var answer = document.querySelector('.answer');
                if(answer != null){ answer.style["display"] = "block"; }
                var message = document.querySelector('.show_massebereich');
                if(message != null){ message.style["display"] = "flex"; }
            }
        },
        error: function(response) {
            $("#.rogressbar").width("0%"); // Handle errors in updating the progress bar
        }
    });
    $.ajax({
        typ: 'GET',
        url: "/webscraper/get_witz", // Request a random message from the server
        success: function(response){
            var show_witz = document.querySelector('.show_witz');
            show_witz.style["display"] = "flex";
            show_witz.innerHTML = response.witz; // Display the random message
        },
        error: function(response){ new_substances_loaded("error"); }, // Handle errors in getting the random message
    });
}

// Function to reset the database
function reset_db(){
    changeBackground(); // Change the background
    var answer = document.querySelector('.answer');
    if(answer != null){ answer.style["display"] = "none"; } // Hide any existing answer
    var message = document.querySelector('.show_massebereich');
    if(message != null){ message.style["display"] = "none"; } // Hide any existing message
    
    // Start showing the progress of the reset
    resetInterval = setInterval(function() {
        show_reset_progress(resetInterval);
    }, 10000); // Update every 10 seconds

    document.querySelector('.progressbar').style["display"] = "flex"; 
    document.querySelector('.progress_message').style["display"] = "flex"; 
    
    $.ajax({
        type: 'GET',
        url: "/webscraper/reset_database", // Request to reset the database
        success: function(response){
            new_substances_loaded("DATABASE WAS COMPLETELY RESTORED"); // Notify the user
            document.querySelector('.show_witz').style["display"] = "none"; // Hide any existing message
        },
        error: function(){
            console.log("error occurred"); // Log any errors
        }
    });
}

// Asynchronous function to load new substances
async function load_new_substances() {
    changeBackground(); // Change the background
    try {
        const response = await $.ajax({
            type: 'GET',
            url: "webscraper/get_new_substances", // Request to get new substances
            timeout: 0
        });

        let message = "New Substances loaded:";
        for (url in response.newSubstances) {
            message += "\n" + response.newSubstances[url]; // Build the message with new substances
        }
        new_substances_loaded(message); // Notify the user
        changeBackground(); // Change the background again
    } catch (error) {
        new_substances_loaded("not possible, try to reload whole db"); // Notify user of error
    }
}

// Function to display a message indicating that new substances are loaded
function new_substances_loaded(message){
    document.querySelector(".container").style["background-color"] = "black"; // Change the background color
    let show_message = document.createElement("div");
    show_message.setAttribute("id","message_bar");
    let skull_video = document.createElement("video");
    skull_video.setAttribute("id","skull_video");
    let message_div = document.createElement("div");
    message_div.setAttribute("id","message_div");
    let skull_source = document.createElement("source");
    skull_source.setAttribute("src","{% static 'videos/roboter.mp4' %}"); // Path to video
    skull_source.setAttribute("type","video/mp4");
    skull_video.append(skull_source);

    show_message.append(skull_video);
    skull_video.autoplay = true;
    skull_video.loop = true;
    let bring_message = document.createElement("div");
    bring_message.innerHTML = message;
    bring_message.setAttribute("id","bring_message");

    show_message.append(skull_video);
    
    let dismis = document.createElement("button");
    dismis.setAttribute("id","dismis_message");
    dismis.innerHTML = "Dismiss Information"; // Button to dismiss the message
    message_div.append(dismis);
    message_div.append(bring_message);
    document.querySelector(".container").append(show_message);
    document.querySelector(".container").append(message_div);
    setTimeout(function() {
        skull_video.style.display = "none"; // Hide video after 1 second
        document.querySelector(".container").style["background-color"] = "transparent"; // Restore background color
        show_message.style["display"] = "none";
    }, 1000);
    dismis.onclick = function() {
        message_div.remove(); // Remove message when button is clicked
        show_message.remove();
    };
}

// Event listener to clear suggestions when focus is lost from the search field
document.querySelector('#search_field').addEventListener("focusout", () => {
    setTimeout(() => $('#suggestionList').empty(), 100);    
})

// Function to display suggestions
function displaySuggestions(suggestions) {
    $('#suggestionList').empty(); // Clear existing suggestions
    for (var suggestion in suggestions) {
        var displayed_text = suggestions[suggestion];
        const sug = $('<div></div>');
        sug.text(displayed_text);
        sug.addClass("suggestionItem");
        sug.on("click", function() {
            var text = $(this).text();
            $('#search_field').val(text); // Set the search field value to the clicked suggestion
            $('#suggestionList').empty(); // Clear suggestions
        });
        $('#suggestionList').append(sug); // Append suggestion to the list
    }
}

// Function to generate a JSON file
function generateJsonFile(){
    try {
        $.ajax({
            type: "GET",
            url: "webscraper/generate", // Request to generate a JSON file
            success: function(response){
                new_substances_loaded("Here is the path to the generated File: " + response); // Notify the user
            }
        });
    } catch (Exception) {
        new_substances_loaded("This action isn't possible"); // Notify user of error
    }
}

// Function to process the JSON input file
function proccesJsonFile(){
    var err = false;
    $.ajax({
        type: "GET",
        url: "webscraper/processJsonInput", // Request to process JSON input
        success: function(response){
            new_substances_loaded(response); // Notify the user
            err = true;
        },
        error: function(response){
            new_substances_loaded("File couldn't be stored"); // Notify user of error
        },
    });
    setTimeout(() => {
        if (err == false) {
            new_substances_loaded("Starting to process the Files in the Directory: files_to_insert\nplease wait until completion\nyou will be informed"); // Notify user of ongoing process
        }
    }, 1500);
}

// Function to get the selected search category
function get_search_category(){
    return $('#select').val();
}

// Function to get the value of the search field
function get_what_is_searched(){
    return $('#search_field').val();
}

const suggestions = [];

// Function to handle search suggestions
function searchSuggestions(){
    var category = get_search_category();
    var what_is_searched = get_what_is_searched();

    // Conditions to request search suggestions
    if ((what_is_searched.length >= 2 && category != "molecular_mass" && category != "smiles") || 
        (what_is_searched.length > 5 && category == "smiles") || 
        category == "category") {

        data = {
            category: category,
            searched: what_is_searched,
        }

        $.ajax({
            type: "GET",
            url: "webscraper/search_suggestion", // Request search suggestions
            data: data,
            success: (response) => {
                let suggestions = response;
                displaySuggestions(suggestions, category); // Display the suggestions
            },
            error: (response) => {
                // Handle errors if needed
            },
        });
    }
}

// Function to delete search results
function delete_search_results() {
    new_substances_loaded("Delete request made"); // Notify the user
    $.ajax({
        type: 'GET',
        url: "webscraper/delete_search_result", // Request to delete search results
        success: function(response) {
            new_substances_loaded("Delete request successful"); // Notify the user
            // Process server response if needed
        },
        error: function(xhr, status, error) {
            console.error("Error deleting search results:", error); // Log any errors
            // Handle errors if needed
        }
    });
}

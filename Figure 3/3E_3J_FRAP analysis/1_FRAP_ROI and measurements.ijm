// Note to update paths at line 18 and 25

// Initialize loop variable
number_of_ROIs = 0;

// Start loop
while (true) {
    // Step 1: Draw an ROI
    run("ROI Manager...");
    waitForUser("Draw ROI, click 'Add' in the ROI Manager.");
    waitForUser("Select the ROI to save.");

    // Increment the number of ROIs
        number_of_ROIs = number_of_ROIs + 1;
        
    // Step 2: Save the ROI
    // Change path ".../ROI_", add file number before "ROI_"
    roiManager("Save", "~/ROI_" + number_of_ROIs + ".roi");

    // Step 3: Measure Stack
    run("Set Measurements...", "area mean standard min integrated redirect=None decimal=3");
    run("Measure Stack...");

    // Step 4: Save Measurements
    // Change path ".../Measurements_", add file number before "ROI_"
    saveAs("Results", "~/Measurements_" + number_of_ROIs + ".csv");

    // Close the "Results" window
    selectWindow("Results");
    run("Close");

    // Ask if the user wants to repeat the process
    yesNo = getBoolean("Do you want to draw another ROI?", "Yes", "No");
    
    // If the user chooses 'No', break the loop
    if (!yesNo) {
        showMessage("Loop Stopped", "Good job!");
        break;
    }
}

// End of the loop

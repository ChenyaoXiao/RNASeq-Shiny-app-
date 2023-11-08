# RNASeq-Shiny-app-
This app was designed to facilitate quick searches for the expression levels of one or more genes in any dataset. 
The user interface is dynamic, simple, and intuitive to ensure a seamless experience. 

## Dataset information
The left side shows the user interface organized by program, model type, model name and tissue type, 
while the right side displays the output figures and data retrieval options. 
The dataset selection drop-down menu changes based on your selection. If you click the checkbox below the radio buttons, 
the menu will display all datasets containing normal tisse such as PBMC or mouse whole blood. 
If you leave it unchecked, all datasets will be visible in the dropdown menu.

If you need to check any details about program, model, and normal tissue message, you can find all the relevant information in the 
browse dataset panel. Additionally, you can search for details using the search box.

## Gene symbol Search

You can enter gene symbols you wish to search for in the text box. If you enter a gene symbol like RAM that is not valid, 
the app will display an error message and ask you to confirm the symbol. The data table will then show all symbols containing RAM, 
which can help suggest the correct symbol you intended to search for. 
You can also enter multiple gene symbols separated by semicolons. After that, you can select a plot type from the options provided, 
which include bar chart, box plot, and violin plot, to visualize the data. When using the data table option, you can choose to display data for TPM and fold change. 

## Excel file download

The app also has a download button for downloading data table results as an excel file according to your choice.

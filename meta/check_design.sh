#!/bin/bash

# Function 1 : Check existing folder

# Set the path to the directory where you want to check for folders
directory="/group/dairy/BovLRC/AgVic/Basecall"

# Set the path to your CSV file
csv_file="/group/dairy/BovLRC/nf-EXPLOR/meta/metadata_LR.csv"

gawk -F, '{print $1}' "$csv_file" > tmp

# Loop through each line in the CSV file
while IFS=, read -r id
do
    # Check if a folder with the name of the ID exists in the directory
    if [ -d "$directory/$id" ]; then
        echo "Folder $id exists"
    else
        echo "Folder $id does not exist"
    fi
done < tmp

rm tmp

# Function 2 : Check POD5 folder

# Get column from user
column_num=4

check_directory() {
    # Open the CSV file
    while IFS=, read -r -a row; do
        # Get the value in the specified column
        value="${row[$column_num - 1]}"
        # Check if the directory exists
        if [ -d "$value" ]; then
            echo "The directory $value exists"
        else
            echo "The directory $value DOES NOT exist"
        fi
    done < "$csv_file"
}

# Check if CSV file exists
if [ ! -f "$csv_file" ]; then
    echo "CSV file not found!"
    exit 1
fi

# Execute the function
check_directory
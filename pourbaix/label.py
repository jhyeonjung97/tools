#!/usr/bin/env python3
import os
import glob
import csv

def create_label_csv():
    """
    Creates a label.csv file with JSON filenames in the current folder listed in the first column,
    followed by five empty columns.
    """
    # Find .json files in the current folder
    json_files = glob.glob("*.json")
    
    # Extract filenames only (without path)
    json_filenames = [os.path.basename(file) for file in json_files]
    
    # Sort filenames (optional)
    json_filenames.sort()
    
    # Create label.csv file
    with open('label.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        if json_filenames:
            # If JSON files exist, write each filename in the first column with five empty columns
            for filename in json_filenames:
                writer.writerow([filename, '', '', '', '', ''])
            print(f"label.csv has been created. Found {len(json_filenames)} JSON files:")
            for filename in json_filenames:
                print(f"  - {filename}")
        else:
            # If no JSON files exist, create an empty row with five empty columns
            writer.writerow(['', '', '', '', '', ''])
            print("No JSON files found in the current folder. An empty label.csv has been created.")

if __name__ == "__main__":
    create_label_csv() 
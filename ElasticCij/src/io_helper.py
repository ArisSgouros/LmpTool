#!/usr/bin/env python3
import os, sys

def replace_string_in_file(file_path, old_string, new_string):
    """
    Replaces all occurrences of old_string with new_string in a file.

    Args:
        file_path (str): The path to the file.
        old_string (str): The string to be replaced.
        new_string (str, int, float): The string to replace with.
    """
    try:
        # Make sure that the input value is a string
        new_string = str(new_string)

        # Read the entire file content into memory
        with open(file_path, 'r') as file:
            file_content = file.read()

        # Perform the replacement
        new_content = file_content.replace(old_string, new_string)

        # Write the modified content back to the file
        with open(file_path, 'w') as file:
            file.write(new_content)

        #print(f"Successfully replaced '{old_string}' with '{new_string}' in '{file_path}'.")

    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

def append_line_to_file(file_path, string_to_append):
    """
    Appends a string to the end of a file.

    Args:
        file_path (str): The path to the file.
        string_to_append (str): The string to be appended.
    """
    try:
        # Open the file in append mode ('a')
        # The 'with' statement ensures the file is automatically closed
        with open(file_path, 'a') as file:
            file.write("%s\n" % string_to_append)
            #print(f"Successfully appended '{string_to_append}' to '{file_path}'.")
    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

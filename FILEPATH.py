# FILEPATH: Untitled-1

# Read the file
with open('C185data_updated.txt', 'r') as file:
    lines = file.readlines()

# Sort the lines based on the third column
sorted_lines = sorted(lines, key=lambda line: line.split(',')[3])

# Write the sorted data back to the file
with open('output2.txt', 'w') as file:
    file.writelines(sorted_lines)

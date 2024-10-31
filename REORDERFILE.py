
# Read the file
with open('output/C185data_1overtone.txt', 'r') as file:
    lines = file.readlines()

# Sort the lines based on the third column
sorted_lines = sorted(lines, key=lambda line: line.split(',')[3])

# Write the sorted data back to the file
with open('output3.txt', 'w') as file:
    file.writelines(sorted_lines)

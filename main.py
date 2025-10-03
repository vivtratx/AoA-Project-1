import csv


# Load all the points into an array of points
points = []

with open('input.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for line in csvFile:
        temp = tuple(line)
        points.append(temp)


# Main Logic goes down here:
print(points)



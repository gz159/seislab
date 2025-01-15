import os

#---------------------------------------------------------------------------------------
def point_in_polygon(x, y, xy):
    n = len(xy) // 2  # Number of points in polygon.
    x_coords = xy[0::2]  # Even indices: x-coordinates.
    y_coords = xy[1::2]  # Odd indices: y-coordinates.

    side = False  # 0 = outside, 1 = inside.

    for i in range(n):
        j = (i + 1) % n

        if (
            ((y_coords[i] <= y < y_coords[j]) or (y_coords[j] <= y < y_coords[i]))
            and (x < (x_coords[j] - x_coords[i]) * (y - y_coords[i]) / (y_coords[j] - y_coords[i]) + x_coords[i])
        ):
            side = not side  # Jump the fence.

    return True if side else False

#---------------------------------------------------------------------------------------
# Read the polygon coordinates from SeleArea.txt
clippath = "SeleArea.txt"
clippathpolygon = []
try:
    with open(clippath, 'r') as file:
        for line in file:
            data = line.split()
            lat, lon = float(data[0]), float(data[1])
            clippathpolygon.extend([lat, lon])
except IOError as e:
    print(f"Can't open select region coordinates file: {clippath}\nError: {e}")
    exit()

#---------------------------------------------------------------------------------------
# Create a new file to keep selected points
deviashearv = "CEICselesta.dat"
if os.path.exists(deviashearv):
    os.remove(deviashearv)

#---------------------------------------------------------------------------------------
# Read all points data from ChinaCEIC.txt and select those within the polygon
inputf = "ChinaCEIC.txt"
print(f"Processing CEIC station file: {inputf}")

try:
    with open(inputf, 'r') as infile, open(deviashearv, 'w') as outfile:
        for line in infile:
            data = line.split()
            if data[0].startswith("#"):
                continue
            else:
                lat, lon = float(data[4]), float(data[5])
                if point_in_polygon(lat, lon, clippathpolygon):
                    outfile.write(line)
except IOError as e:
    print(f"Can't open CEIC station file: {inputf}\nError: {e}")
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

    return 1 if side else 0


# Example usage:
polygon = [1, 1, 3, 5, 6, 2, 9, 6, 10, 0, 4, 2, 5, -2]
print("(2, 4):", point_in_polygon(2, 4, polygon))  # Expected: 1 (inside)
print("(3, 0.5):", point_in_polygon(3, 0.5, polygon))  # Expected: 0 (outside)
print("(5, 1.5):", point_in_polygon(5, 1.5, polygon))  # Expected: 0 (outside)




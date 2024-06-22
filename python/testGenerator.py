import numpy as np
import csv
from mpmath import mp, cos, sin

mp.dps = 150

num_points = 100
random_points = np.random.uniform(-1.0, 1.0, num_points)

results_cos = [(x, cos(x)) for x in random_points]
with open('../test-data/cos_points.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['x', 'cos(x)'])
    for x, cos_x in results_cos:
        csvwriter.writerow([str(x), str(cos_x)])

results_sin = [(x, sin(x)) for x in random_points]
with open('../test-data/sin_points.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['x', 'sin(x)'])
    for x, sin_x in results_sin:
        csvwriter.writerow([str(x), str(sin_x)])

print("Результаты сохранены в файлы")
